#include <gss/clique.hh>
#include <gss/configuration.hh>
#include <gss/homomorphism.hh>
#include <gss/innards/graph_traits.hh>
#include <gss/innards/homomorphism_domain.hh>
#include <gss/innards/homomorphism_model.hh>
#include <gss/innards/homomorphism_proofs.hh>
#include <gss/innards/homomorphism_searcher.hh>
#include <gss/innards/homomorphism_traits.hh>
#include <gss/innards/proof.hh>
#include <gss/innards/solve_state.hh>
#include <gss/innards/thread_utils.hh>

#include <algorithm>
#include <atomic>
#include <barrier>
#include <condition_variable>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <thread>
#include <unordered_set>
#include <utility>

#include <gss/utils/hashing_utils.hh>

using namespace gss;
using namespace gss::innards;

using std::atomic;
using std::barrier;
using std::function;
using std::make_optional;
using std::make_shared;
using std::make_unique;
using std::map;
using std::move;
using std::mutex;
using std::optional;
using std::pair;
using std::shared_ptr;
using std::size_t;
using std::sort;
using std::string;
using std::thread;
using std::to_string;
using std::unique_lock;
using std::unique_ptr;
using std::unordered_set;
using std::vector;

using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::chrono::steady_clock;
using std::chrono::operator""ms;

namespace
{
    struct HomomorphismSolver
    {
        using Domains = vector<HomomorphismDomain>;

        SolveState & state;
        const HomomorphismModel & model;
        const HomomorphismParams & params;
        const std::shared_ptr<Proof> proof;

        HomomorphismSolver(SolveState & s, const HomomorphismParams & p,
            const std::shared_ptr<Proof> & r) :
            state(s),
            model(*s.model),
            params(p),
            proof(r)
        {
        }
    };

    struct SequentialSolver : HomomorphismSolver
    {
        using HomomorphismSolver::HomomorphismSolver;

        auto solve() -> HomomorphismResult
        {
            HomomorphismResult result;

            // domains (carried in the shared solve state)
            Domains & domains = state.domains;
            domains = Domains(model.pattern_size, HomomorphismDomain{model.target_size});
            if (! model.initialise_domains(domains)) {
                result.complete = true;
                model.add_extra_stats(result.extra_stats);
                return result;
            }

            // assignments
            HomomorphismAssignments assignments;
            assignments.values.reserve(model.pattern_size);

            // start search timer
            auto search_start_time = steady_clock::now();

            // do the search
            bool done = false;
            unsigned number_of_restarts = 0;

            HomomorphismSearcher searcher(
                model, params, [](const HomomorphismAssignments &) -> bool { return true; }, proof, state.watches);

            while (! done) {
                ++number_of_restarts;

                // start watching new nogoods
                done = searcher.watches.apply_new_nogoods(
                    [&](const HomomorphismAssignment & assignment) {
                        for (auto & d : domains)
                            if (d.v == assignment.pattern_vertex) {
                                d.values.reset(assignment.target_vertex);
                                d.count = d.values.count();
                                done = done || (0 == d.count);
                                break;
                            }
                    });

                if (done) {
                    result.complete = true;
                    break;
                }

                searcher.watches.clear_new_nogoods();

                ++result.propagations;
                if (searcher.propagate(true, domains, assignments)) {
                    auto assignments_copy = assignments;

                    switch (searcher.restarting_search(assignments_copy, domains, result.nodes, result.propagations,
                        result.solution_count, 0, *params.restarts_schedule)) {
                    case SearchResult::Satisfiable:
                        searcher.save_result(assignments_copy, result);
                        // when counting, reaching here means the enumerate callback asked
                        // us to stop (e.g. a solution limit), so the search is not complete
                        result.complete = ! params.count_solutions;
                        done = true;
                        break;

                    case SearchResult::SatisfiableButKeepGoing:
                        result.complete = true;
                        done = true;
                        break;

                    case SearchResult::Unsatisfiable:
                        result.complete = true;
                        done = true;
                        break;

                    case SearchResult::Aborted:
                        done = true;
                        break;

                    case SearchResult::Restart:
                        break;
                    }
                }
                else {
                    if (proof)
                        proof->root_propagation_failed();
                    result.complete = true;
                    done = true;
                }

                params.restarts_schedule->did_a_restart();
            }

            if (params.restarts_schedule->might_restart())
                result.extra_stats.emplace_back("restarts = " + to_string(number_of_restarts));

            result.extra_stats.emplace_back("shape_graphs = " + to_string(model.max_graphs));

            result.extra_stats.emplace_back("search_time = " + to_string(duration_cast<milliseconds>(steady_clock::now() - search_start_time).count()));

            if (might_have_watches(params)) {
                result.extra_stats.emplace_back("nogoods_size = " + to_string(searcher.watches.nogoods.size()));

                map<int, int> nogoods_lengths;
                for (auto & n : searcher.watches.nogoods)
                    nogoods_lengths[n.literals.size()]++;

                string nogoods_lengths_str;
                for (auto & n : nogoods_lengths) {
                    nogoods_lengths_str += " ";
                    nogoods_lengths_str += to_string(n.first) + ":" + to_string(n.second);
                }
                result.extra_stats.emplace_back("nogoods_lengths =" + nogoods_lengths_str);
            }

            model.add_extra_stats(result.extra_stats);
            return result;
        }
    };

    struct VertexToVertexMappingHash
    {
        auto operator()(const VertexToVertexMapping & v) const -> size_t
        {
            size_t result{0};
            for (auto & [p, t] : v) {
                hash_combine(result, p);
                hash_combine(result, t);
            }
            return result;
        }
    };

    struct ThreadedSolver : HomomorphismSolver
    {
        unsigned n_threads;

        ThreadedSolver(SolveState & s, const HomomorphismParams & p,
            const std::shared_ptr<Proof> & r, unsigned t) :
            HomomorphismSolver(s, p, r),
            n_threads(t)
        {
        }

        auto solve() -> HomomorphismResult
        {
            mutex common_result_mutex;
            HomomorphismResult common_result;
            string by_thread_nodes, by_thread_propagations;

            // domains (the canonical root domains, carried in the shared solve state;
            // each thread takes its own working copy below)
            Domains & common_domains = state.domains;
            common_domains = Domains(model.pattern_size, HomomorphismDomain{model.target_size});
            if (! model.initialise_domains(common_domains)) {
                common_result.complete = true;
                return common_result;
            }

            // start search timer
            auto search_start_time = steady_clock::now();

            vector<thread> threads;
            threads.reserve(n_threads);

            vector<unique_ptr<HomomorphismSearcher>> searchers{n_threads};

            // each thread keeps its own nogood store (the threaded search is the
            // terminal, unbounded step, so these are not the carried state.watches)
            vector<Watches<HomomorphismAssignment, HomomorphismAssignmentWatchTable>> thread_watches{n_threads};

            barrier wait_for_new_nogoods_barrier{n_threads}, synced_nogoods_barrier{n_threads};
            atomic<bool> restart_synchroniser{false};

            mutex duplicate_filter_set_mutex;
            unordered_set<VertexToVertexMapping, VertexToVertexMappingHash> duplicate_filter_set;

            function<auto(unsigned)->void> work_function =
                [&searchers, &common_domains, &thread_watches, &threads, &work_function,
                    &model = this->model, &params = this->params, proof = this->proof, n_threads = this->n_threads,
                    &common_result, &common_result_mutex, &by_thread_nodes, &by_thread_propagations,
                    &wait_for_new_nogoods_barrier, &synced_nogoods_barrier, &restart_synchroniser,
                    &duplicate_filter_set, &duplicate_filter_set_mutex](unsigned t) -> void {
                // do the search
                HomomorphismResult thread_result;

                bool just_the_first_thread = (0 == t) && params.delay_thread_creation;

                searchers[t] = make_unique<HomomorphismSearcher>(
                    model, params, [&](const HomomorphismAssignments & a) -> bool {
                        VertexToVertexMapping v;
                        searchers[t]->expand_to_full_result(a, v);
                        unique_lock<mutex> lock{duplicate_filter_set_mutex};
                        return duplicate_filter_set.insert(v).second;
                    },
                    proof, thread_watches[t]);
                if (0 != t)
                    searchers[t]->set_seed(t);

                unsigned number_of_restarts = 0;

                Domains domains = common_domains;

                HomomorphismAssignments thread_assignments;
                thread_assignments.values.reserve(model.pattern_size);

                // each thread needs its own restarts schedule
                unique_ptr<RestartsSchedule> thread_restarts_schedule;
                if (0 == t || ! params.triggered_restarts)
                    thread_restarts_schedule.reset(params.restarts_schedule->clone());
                else
                    thread_restarts_schedule = make_unique<SyncedRestartSchedule>(restart_synchroniser);

                while (true) {
                    ++number_of_restarts;

                    if (! just_the_first_thread) {
                        wait_for_new_nogoods_barrier.arrive_and_wait();

                        for (unsigned u = 0; u < n_threads; ++u)
                            if (t != u)
                                searchers[t]->watches.gather_nogoods_from(searchers[u]->watches);

                        // start watching new nogoods
                        if (searchers[t]->watches.apply_new_nogoods(
                                [&](const HomomorphismAssignment & assignment) {
                                    for (auto & d : domains)
                                        if (d.v == assignment.pattern_vertex) {
                                            d.values.reset(assignment.target_vertex);
                                            d.count = d.values.count();
                                            break;
                                        }
                                }))
                            break;

                        if (0 == t) {
                            restart_synchroniser.store(false);
                            duplicate_filter_set.clear();
                        }

                        synced_nogoods_barrier.arrive_and_wait();

                        searchers[t]->watches.clear_new_nogoods();
                    }

                    ++thread_result.propagations;
                    if (searchers[t]->propagate(true, domains, thread_assignments)) {
                        auto assignments_copy = thread_assignments;

                        switch (searchers[t]->restarting_search(assignments_copy, domains, thread_result.nodes, thread_result.propagations,
                            thread_result.solution_count, 0, *thread_restarts_schedule)) {
                        case SearchResult::Satisfiable:
                            searchers[t]->save_result(assignments_copy, thread_result);
                            thread_result.complete = true;
                            params.timeout->trigger_early_abort();
                            searchers[t]->watches.post_nogood(Nogood<HomomorphismAssignment>{});
                            break;

                        case SearchResult::SatisfiableButKeepGoing:
                            thread_result.complete = true;
                            params.timeout->trigger_early_abort();
                            searchers[t]->watches.post_nogood(Nogood<HomomorphismAssignment>{});
                            break;

                        case SearchResult::Unsatisfiable:
                            thread_result.complete = true;
                            params.timeout->trigger_early_abort();
                            searchers[t]->watches.post_nogood(Nogood<HomomorphismAssignment>{});
                            break;

                        case SearchResult::Aborted:
                            searchers[t]->watches.post_nogood(Nogood<HomomorphismAssignment>{});
                            break;

                        case SearchResult::Restart:
                            break;
                        }
                    }
                    else {
                        thread_result.complete = true;
                        searchers[t]->watches.post_nogood(Nogood<HomomorphismAssignment>{});
                        params.timeout->trigger_early_abort();
                    }

                    if (0 == t)
                        restart_synchroniser.store(true);
                    thread_restarts_schedule->did_a_restart();

                    if (params.delay_thread_creation && just_the_first_thread) {
                        if (! thread_result.complete) {
                            just_the_first_thread = false;
                            for (unsigned u = 1; u < n_threads; ++u)
                                threads.emplace_back([&, u]() { work_function(u); });
                        }
                        else
                            break;
                    }
                }

                if (params.delay_thread_creation && 0 == t)
                    for (auto & th : threads)
                        th.join();

                unique_lock<mutex> lock{common_result_mutex};
                if (! thread_result.mapping.empty())
                    common_result.mapping = move(thread_result.mapping);
                common_result.nodes += thread_result.nodes;
                common_result.propagations += thread_result.propagations;
                common_result.solution_count += thread_result.solution_count;
                common_result.complete = common_result.complete || thread_result.complete;
                for (auto & x : thread_result.extra_stats)
                    common_result.extra_stats.push_back("t" + to_string(t) + "_" + x);

                by_thread_nodes.append(" " + to_string(thread_result.nodes));
                by_thread_propagations.append(" " + to_string(thread_result.propagations));
            };

            if (params.delay_thread_creation)
                work_function(0);
            else {
                for (unsigned u = 0; u < n_threads; ++u)
                    threads.emplace_back([&, u]() { work_function(u); });

                for (auto & th : threads)
                    th.join();
            }

            common_result.extra_stats.emplace_back("by_thread_nodes =" + by_thread_nodes);
            common_result.extra_stats.emplace_back("by_thread_propagations =" + by_thread_propagations);
            common_result.extra_stats.emplace_back("search_time = " + to_string(duration_cast<milliseconds>(steady_clock::now() - search_start_time).count()));

            return common_result;
        }
    };

    // --- Solve pipeline -------------------------------------------------------------
    //
    // The top-level solve is expressed as a sequence of steps over a shared context.
    // Each step may transform the context and may conclude the problem (emitting the
    // proof's conclusion); the driver runs them in registration order until one
    // concludes, with the search as the terminal step. This is the skeleton of the
    // unified pipeline described in dev_docs/preprocessor-refactor.md; later phases fold
    // the model, domains and watches into a single SolveState and turn the preprocessing
    // filters and supplemental-graph builders into steps too.

    enum class StepOutcome
    {
        Continue, // carry on to the next step
        Concluded // the problem is decided; stop the pipeline
    };

    struct SolveContext
    {
        const InputGraph & pattern;
        const InputGraph & target;
        const HomomorphismParams & params;
        const shared_ptr<Proof> & proof;
        HomomorphismProofs * hom_proofs = nullptr;
        SolveState state;
        HomomorphismResult result;
    };

    struct SolveStep
    {
        virtual ~SolveStep() = default;
        virtual auto run(SolveContext & ctx) -> StepOutcome = 0;
    };

    // Emit the OPB model through the solver-proofs layer (variables, (local-)injectivity,
    // adjacency / induced non-edges, the preserved set, then finalise + loop-fix). Runs
    // first so the later steps' proof conclusions can reference the model. No-op when
    // proof logging is off.
    struct EmitProofModelStep : SolveStep
    {
        auto run(SolveContext & ctx) -> StepOutcome override
        {
            if (ctx.hom_proofs)
                ctx.hom_proofs->emit_model(ctx.pattern, ctx.target, ctx.params);
            return StepOutcome::Continue;
        }
    };

    // If we are finding a non-shrinking mapping and the pattern has more vertices than
    // the target, there is trivially no solution.
    struct PatternBiggerThanTargetStep : SolveStep
    {
        auto run(SolveContext & ctx) -> StepOutcome override
        {
            if (! (is_nonshrinking(ctx.params) && (ctx.pattern.size() > ctx.target.size())))
                return StepOutcome::Continue;

            if (ctx.proof) {
                ctx.proof->failure_due_to_pattern_bigger_than_target();
                if (ctx.params.count_solutions)
                    ctx.proof->finish_enumeration_proof(0, true);
                else
                    ctx.proof->finish_unsat_proof();
            }

            ctx.result = HomomorphismResult{};
            return StepOutcome::Concluded;
        }
    };

    // If the target has a loop and we only need a single non-injective mapping, map
    // every pattern vertex onto the loop.
    struct TargetLoopShortcutStep : SolveStep
    {
        auto run(SolveContext & ctx) -> StepOutcome override
        {
            const auto & params = ctx.params;
            const auto & pattern = ctx.pattern;
            const auto & target = ctx.target;

            if (! ((params.injectivity == Injectivity::NonInjective) && ! pattern.has_vertex_labels() && ! pattern.has_edge_labels() && target.loopy() && ! params.count_solutions && ! params.enumerate_callback))
                return StepOutcome::Continue;

            HomomorphismResult result;
            result.extra_stats.emplace_back("used_loops_property = true");
            result.complete = true;
            int loop = -1;
            for (int t = 0; t < target.size(); ++t)
                if (target.adjacent(t, t)) {
                    loop = t;
                    break;
                }

            for (int n = 0; n < pattern.size(); ++n)
                result.mapping.emplace(n, loop);

            ++result.solution_count;

            ctx.result = move(result);
            return StepOutcome::Concluded;
        }
    };

    // If the pattern is a clique, solve it with the clique algorithm instead.
    struct CliqueShortcutStep : SolveStep
    {
        auto run(SolveContext & ctx) -> StepOutcome override
        {
            const auto & params = ctx.params;
            const auto & pattern = ctx.pattern;
            const auto & target = ctx.target;

            if (! (can_use_clique(params) && is_simple_clique(pattern)))
                return StepOutcome::Continue;

            CliqueParams clique_params;
            clique_params.timeout = params.timeout;
            clique_params.start_time = params.start_time;
            clique_params.decide = make_optional(pattern.size());
            clique_params.restarts_schedule = make_unique<NoRestartsSchedule>();
            auto clique_result = solve_clique_problem(target, clique_params);

            // now translate the result back into what we expect
            HomomorphismResult result;
            int v = 0;
            for (auto & m : clique_result.clique) {
                result.mapping.emplace(v++, m);
                // the clique solver can find a bigger clique than we ask for
                if (v >= pattern.size())
                    break;
            }
            result.nodes = clique_result.nodes;
            result.extra_stats = move(clique_result.extra_stats);
            result.extra_stats.emplace_back("used_clique_solver = true");
            result.complete = clique_result.complete;

            ctx.result = move(result);
            return StepOutcome::Concluded;
        }
    };

    // Build the constraint model and run the search. The terminal step: it always
    // concludes (with a solution, unsat, or an aborted/partial result).
    struct MainSolveStep : SolveStep
    {
        auto run(SolveContext & ctx) -> StepOutcome override
        {
            const auto & pattern = ctx.pattern;
            const auto & target = ctx.target;
            const auto & params = ctx.params;
            const auto & proof = ctx.proof;

            // Build the model now -- after the cheap shortcut steps have had their
            // chance to conclude -- into the carried state, where the search steps
            // (and, later, staged builder / filter steps) read it.
            ctx.state.model = make_unique<HomomorphismModel>(target, pattern, params, proof, ctx.hom_proofs);
            auto & model = *ctx.state.model;

            if (! model.prepare()) {
                HomomorphismResult result;
                result.extra_stats.emplace_back("model_consistent = false");
                result.complete = true;
                if (proof) {
                    if (params.count_solutions)
                        proof->finish_enumeration_proof(0, true);
                    else
                        proof->finish_unsat_proof();
                }
                ctx.result = move(result);
                return StepOutcome::Concluded;
            }

            HomomorphismResult result;
            if (1 == params.n_threads) {
                SequentialSolver solver(ctx.state, params, proof);
                result = solver.solve();
            }
            else {
                if (! params.restarts_schedule->might_restart())
                    throw UnsupportedConfiguration{"Threaded search requires restarts"};

                unsigned n_threads = how_many_threads(params.n_threads);
                ThreadedSolver solver(ctx.state, params, proof, n_threads);
                result = solver.solve();
            }

            if (proof) {
                if (params.count_solutions)
                    // counting / enumeration: a complete search yields ENUMERATION_COMPLETE,
                    // otherwise (timeout or solution limit) ENUMERATION_PARTIAL
                    proof->finish_enumeration_proof(result.solution_count, result.complete);
                else if (result.complete && result.mapping.empty())
                    proof->finish_unsat_proof();
                else if (! result.mapping.empty())
                    proof->finish_sat_proof();
                else
                    proof->finish_unknown_proof();
            }

            ctx.result = move(result);
            return StepOutcome::Concluded;
        }
    };
}

auto gss::solve_homomorphism_problem(
    const InputGraph & pattern,
    const InputGraph & target,
    const HomomorphismParams & params) -> HomomorphismResult
{
    // start by setting up proof logging, if necessary
    shared_ptr<Proof> proof;
    if (params.proof_options) {
        // proof logging is currently incompatible with a whole load of "extra" features,
        // but can be adapted to support most of them
        if (1 != params.n_threads)
            throw UnsupportedConfiguration{"Proof logging cannot yet be used with threads"};
        if (params.clique_detection)
            throw UnsupportedConfiguration{"Proof logging cannot yet be used with clique detection, use --no-clique-detection"};
        if (! params.pattern_less_constraints.empty() || ! params.target_occur_less_constraints.empty())
            throw UnsupportedConfiguration{"Proof logging cannot yet be used with less-constraints"};
        if (pattern.has_vertex_labels() || pattern.has_edge_labels())
            throw UnsupportedConfiguration{"Proof logging cannot yet be used on labelled graphs"};
        if (params.count_solutions && params.restarts_schedule && params.restarts_schedule->might_restart())
            throw UnsupportedConfiguration{"Proof logging cannot yet be used when counting with restarts, use --restarts none"};
        proof = make_shared<Proof>(*params.proof_options);
    }

    // The solver-proofs middle layer (owns vertex naming + the homomorphism-specific
    // derivations). Created here so it spans the whole solve: the OPB model emission
    // step below and the model's later supplemental / filter proofs share this one
    // instance. Null when proof logging is off.
    unique_ptr<HomomorphismProofs> hom_proofs;
    if (proof)
        hom_proofs = make_unique<HomomorphismProofs>(proof, pattern, target);

    // The solve runs as a pipeline of steps over a shared context, in registration
    // order, stopping at the first step that concludes the problem; the search is the
    // terminal step (see dev_docs/preprocessor-refactor.md).
    SolveContext ctx{pattern, target, params, proof, hom_proofs.get(), {}, {}};

    vector<unique_ptr<SolveStep>> steps;
    steps.push_back(make_unique<EmitProofModelStep>());          // emit the OPB model
    steps.push_back(make_unique<PatternBiggerThanTargetStep>()); // trivial size refutation
    steps.push_back(make_unique<TargetLoopShortcutStep>());      // non-injective target-loop shortcut
    steps.push_back(make_unique<CliqueShortcutStep>());          // clique-pattern reduction
    steps.push_back(make_unique<MainSolveStep>());               // build model + search (terminal)

    for (auto & step : steps)
        if (step->run(ctx) == StepOutcome::Concluded)
            break;

    return move(ctx.result);
}
