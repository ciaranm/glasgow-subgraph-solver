/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "homomorphism.hh"
#include "clique.hh"
#include "configuration.hh"
#include "graph_traits.hh"
#include "homomorphism_domain.hh"
#include "homomorphism_model.hh"
#include "homomorphism_searcher.hh"
#include "homomorphism_traits.hh"
#include "thread_utils.hh"
#include "proof.hh"

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <thread>
#include <unordered_set>
#include <utility>

#include <boost/thread/barrier.hpp>
#include <boost/functional/hash.hpp>

using std::atomic;
using std::function;
using std::make_optional;
using std::make_unique;
using std::map;
using std::move;
using std::mutex;
using std::optional;
using std::pair;
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

using boost::barrier;
using boost::hash_combine;

namespace
{
    struct HomomorphismSolver
    {
        using Domains = vector<HomomorphismDomain>;

        const HomomorphismModel & model;
        const HomomorphismParams & params;

        HomomorphismSolver(const HomomorphismModel & m, const HomomorphismParams & p) :
            model(m),
            params(p)
        {
        }
    };

    struct SequentialSolver :
        HomomorphismSolver
    {
        using HomomorphismSolver::HomomorphismSolver;

        auto solve() -> HomomorphismResult
        {
            HomomorphismResult result;

            // domains
            Domains domains(model.pattern_size, HomomorphismDomain{ model.target_size });
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

            HomomorphismSearcher searcher(model, params, [] (const HomomorphismAssignments &) -> bool { return true; });

            while (! done) {
                ++number_of_restarts;

                // start watching new nogoods
                done = searcher.watches.apply_new_nogoods(
                        [&] (const HomomorphismAssignment & assignment) {
                            for (auto & d : domains)
                                if (d.v == assignment.pattern_vertex) {
                                    d.values.reset(assignment.target_vertex);
                                    d.count = d.values.count();
                                    done = done || (0 == d.count);
                                    break;
                                }
                        });

                if (done)
                    break;

                searcher.watches.clear_new_nogoods();

                ++result.propagations;
                if (searcher.propagate(true, domains, assignments, params.propagate_using_lackey != PropagateUsingLackey::Never)) {
                    auto assignments_copy = assignments;

                    switch (searcher.restarting_search(assignments_copy, domains, result.nodes, result.propagations,
                                result.solution_count, 0, *params.restarts_schedule)) {
                        case SearchResult::Satisfiable:
                            searcher.save_result(assignments_copy, result);
                            result.complete = true;
                            done = true;
                            break;

                        case SearchResult::SatisfiableButKeepGoing:
                            result.complete = true;
                            done = true;
                            break;

                        case SearchResult::UnsatisfiableAndBackjumpUsingLackey:
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
                    if (params.proof)
                        params.proof->root_propagation_failed();
                    result.complete = true;
                    done = true;
                }

                params.restarts_schedule->did_a_restart();
            }

            if (params.restarts_schedule->might_restart())
                result.extra_stats.emplace_back("restarts = " + to_string(number_of_restarts));

            result.extra_stats.emplace_back("shape_graphs = " + to_string(model.max_graphs));

            result.extra_stats.emplace_back("search_time = " + to_string(
                        duration_cast<milliseconds>(steady_clock::now() - search_start_time).count()));

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
        auto operator() (const VertexToVertexMapping & v) const -> size_t
        {
            size_t result{ 0 };
            for (auto & [ p, t ] : v) {
                hash_combine(result, p);
                hash_combine(result, t);
            }
            return result;
        }
    };

    struct ThreadedSolver : HomomorphismSolver
    {
        unsigned n_threads;

        ThreadedSolver(const HomomorphismModel & m, const HomomorphismParams & p, unsigned t) :
            HomomorphismSolver(m, p),
            n_threads(t)
        {
        }

        auto solve() -> HomomorphismResult
        {
            mutex common_result_mutex;
            HomomorphismResult common_result;
            string by_thread_nodes, by_thread_propagations;

            // domains
            Domains common_domains(model.pattern_size, HomomorphismDomain{ model.target_size });
            if (! model.initialise_domains(common_domains)) {
                common_result.complete = true;
                return common_result;
            }

            // start search timer
            auto search_start_time = steady_clock::now();

            vector<thread> threads;
            threads.reserve(n_threads);

            vector<unique_ptr<HomomorphismSearcher> > searchers{ n_threads };

            barrier wait_for_new_nogoods_barrier{ n_threads }, synced_nogoods_barrier{ n_threads };
            atomic<bool> restart_synchroniser{ false };

            mutex duplicate_filter_set_mutex;
            unordered_set<VertexToVertexMapping, VertexToVertexMappingHash> duplicate_filter_set;

            function<auto (unsigned) -> void> work_function = [&searchers, &common_domains, &threads, &work_function,
                        &model = this->model, &params = this->params, n_threads = this->n_threads,
                        &common_result, &common_result_mutex, &by_thread_nodes, &by_thread_propagations,
                        &wait_for_new_nogoods_barrier, &synced_nogoods_barrier, &restart_synchroniser,
                        &duplicate_filter_set, &duplicate_filter_set_mutex] (unsigned t) -> void
            {
                // do the search
                HomomorphismResult thread_result;

                bool just_the_first_thread = (0 == t) && params.delay_thread_creation;

                searchers[t] = make_unique<HomomorphismSearcher>(model, params, [&] (const HomomorphismAssignments & a) -> bool {
                        VertexToVertexMapping v;
                        searchers[t]->expand_to_full_result(a, v);
                        unique_lock<mutex> lock{ duplicate_filter_set_mutex };
                        return duplicate_filter_set.insert(v).second;
                        });
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
                        wait_for_new_nogoods_barrier.wait();

                        for (unsigned u = 0 ; u < n_threads ; ++u)
                            if (t != u)
                                searchers[t]->watches.gather_nogoods_from(searchers[u]->watches);

                        // start watching new nogoods
                        if (searchers[t]->watches.apply_new_nogoods(
                                [&] (const HomomorphismAssignment & assignment) {
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

                        synced_nogoods_barrier.wait();

                        searchers[t]->watches.clear_new_nogoods();
                    }

                    ++thread_result.propagations;
                    if (searchers[t]->propagate(true, domains, thread_assignments, params.propagate_using_lackey != PropagateUsingLackey::Never)) {
                        auto assignments_copy = thread_assignments;

                        switch (searchers[t]->restarting_search(assignments_copy, domains, thread_result.nodes, thread_result.propagations,
                                    thread_result.solution_count, 0, *thread_restarts_schedule)) {
                            case SearchResult::Satisfiable:
                                searchers[t]->save_result(assignments_copy, thread_result);
                                thread_result.complete = true;
                                params.timeout->trigger_early_abort();
                                searchers[t]->watches.post_nogood(Nogood<HomomorphismAssignment>{ });
                                break;

                            case SearchResult::SatisfiableButKeepGoing:
                                thread_result.complete = true;
                                params.timeout->trigger_early_abort();
                                searchers[t]->watches.post_nogood(Nogood<HomomorphismAssignment>{ });
                                break;

                            case SearchResult::Unsatisfiable:
                            case SearchResult::UnsatisfiableAndBackjumpUsingLackey:
                                thread_result.complete = true;
                                params.timeout->trigger_early_abort();
                                searchers[t]->watches.post_nogood(Nogood<HomomorphismAssignment>{ });
                                break;

                            case SearchResult::Aborted:
                                searchers[t]->watches.post_nogood(Nogood<HomomorphismAssignment>{ });
                                break;

                            case SearchResult::Restart:
                                break;
                        }
                    }
                    else {
                        thread_result.complete = true;
                        searchers[t]->watches.post_nogood(Nogood<HomomorphismAssignment>{ });
                        params.timeout->trigger_early_abort();
                    }

                    if (0 == t)
                        restart_synchroniser.store(true);
                    thread_restarts_schedule->did_a_restart();

                    if (params.delay_thread_creation && just_the_first_thread) {
                        if (! thread_result.complete) {
                            just_the_first_thread = false;
                            for (unsigned u = 1 ; u < n_threads ; ++u)
                                threads.emplace_back([&, u] () { work_function(u); });
                        }
                        else
                            break;
                    }
                }

                if (params.delay_thread_creation && 0 == t)
                    for (auto & th : threads)
                        th.join();

                unique_lock<mutex> lock{ common_result_mutex };
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
                for (unsigned u = 0 ; u < n_threads ; ++u)
                    threads.emplace_back([&, u] () { work_function(u); });

                for (auto & th : threads)
                    th.join();
            }

            common_result.extra_stats.emplace_back("by_thread_nodes =" + by_thread_nodes);
            common_result.extra_stats.emplace_back("by_thread_propagations =" + by_thread_propagations);
            common_result.extra_stats.emplace_back("search_time = " + to_string(
                        duration_cast<milliseconds>(steady_clock::now() - search_start_time).count()));

            return common_result;
        }
    };
}

auto solve_homomorphism_problem(
        const InputGraph & pattern,
        const InputGraph & target,
        const HomomorphismParams & params) -> HomomorphismResult
{
    // start by setting up proof logging, if necessary
    if (params.proof) {
        // proof logging is currently incompatible with a whole load of "extra" features,
        // but can be adapted to support most of them
        if (1 != params.n_threads)
            throw UnsupportedConfiguration{ "Proof logging cannot yet be used with threads" };
        if (params.clique_detection)
            throw UnsupportedConfiguration{ "Proof logging cannot yet be used with clique detection" };
        if (params.lackey)
            throw UnsupportedConfiguration{ "Proof logging cannot yet be used with a lackey" };
        if (! params.pattern_less_constraints.empty() || ! params.target_occur_less_constraints.empty())
            throw UnsupportedConfiguration{ "Proof logging cannot yet be used with less-constraints" };
        if (params.injectivity != Injectivity::Injective && params.injectivity != Injectivity::NonInjective)
            throw UnsupportedConfiguration{ "Proof logging can currently only be used with injectivity or non-injectivity" };
        if (pattern.has_vertex_labels() || pattern.has_edge_labels())
            throw UnsupportedConfiguration{ "Proof logging cannot yet be used on labelled graphs" };

        // set up our model file, with a set of OPB variables for each CP variable
        for (int n = 0 ; n < pattern.size() ; ++n) {
            params.proof->create_cp_variable(n, target.size(),
                    [&] (int v) { return pattern.vertex_name(v); },
                    [&] (int v) { return target.vertex_name(v); });
        }

        // generate constraints for injectivity
        if (params.injectivity == Injectivity::Injective)
            params.proof->create_injectivity_constraints(pattern.size(), target.size());

        // generate edge constraints, and also handle loops here
        for (int p = 0 ; p < pattern.size() ; ++p) {
            for (int t = 0 ; t < target.size() ; ++t) {
                if (pattern.adjacent(p, p) && ! target.adjacent(t, t))
                    params.proof->create_forbidden_assignment_constraint(p, t);
                else if (params.induced && ! pattern.adjacent(p, p) && target.adjacent(t, t))
                    params.proof->create_forbidden_assignment_constraint(p, t);

                // it's simpler to always have the adjacency constraints, even
                // if the assignment is forbidden
                params.proof->start_adjacency_constraints_for(p, t);

                // if p can be mapped to t, then each neighbour of p...
                for (int q = 0 ; q < pattern.size() ; ++q)
                    if (q != p && pattern.adjacent(p, q)) {
                        // ... must be mapped to a neighbour of t
                        vector<int> permitted;
                        for (int u = 0 ; u < target.size() ; ++u)
                            if (t != u && target.adjacent(t, u))
                                permitted.push_back(u);
                        params.proof->create_adjacency_constraint(p, q, t, permitted, false);
                    }

                // same for non-adjacency for induced
                if (params.induced) {
                    for (int q = 0 ; q < pattern.size() ; ++q)
                        if (q != p && ! pattern.adjacent(p, q)) {
                            // ... must be mapped to a neighbour of t
                            vector<int> permitted;
                            for (int u = 0 ; u < target.size() ; ++u)
                                if (t != u && ! target.adjacent(t, u))
                                    permitted.push_back(u);
                            params.proof->create_adjacency_constraint(p, q, t, permitted, true);
                        }
                }
            }
        }

        // output the model file
        params.proof->finalise_model();
    }

    // first sanity check: if we're finding an injective mapping, and there
    // aren't enough vertices, fail immediately.
    if (is_nonshrinking(params) && (pattern.size() > target.size())) {
        if (params.proof) {
            params.proof->failure_due_to_pattern_bigger_than_target();
            params.proof->finish_unsat_proof();
        }

        return HomomorphismResult{ };
    }

    // does the target have loops, and are we looking for a single non-injective mapping?
    if ((params.injectivity == Injectivity::NonInjective) && ! pattern.has_vertex_labels()
           && ! pattern.has_edge_labels() && target.loopy() && ! params.count_solutions
           && ! params.enumerate_callback) {
        HomomorphismResult result;
        result.extra_stats.emplace_back("used_loops_property = true");
        result.complete = true;
        int loop = -1;
        for (int t = 0 ; t < target.size() ; ++t)
            if (target.adjacent(t, t)) {
                loop = t;
                break;
            }

        for (int n = 0 ; n < pattern.size() ; ++n)
            result.mapping.emplace(n, loop);

        return result;
    }

    // is the pattern a clique? if so, use a clique algorithm instead
    if (can_use_clique(params) && is_simple_clique(pattern)) {
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

        return result;
    }
    else {
        // just solve the problem
        HomomorphismModel model(target, pattern, params);

        if (! model.prepare()) {
            HomomorphismResult result;
            result.extra_stats.emplace_back("model_consistent = false");
            result.complete = true;
            if (params.proof)
                params.proof->finish_unsat_proof();
            return result;
        }

        HomomorphismResult result;
        if (1 == params.n_threads) {
            SequentialSolver solver(model, params);
            result = solver.solve();
        }
        else {
            if (! params.restarts_schedule->might_restart())
                throw UnsupportedConfiguration{ "Threaded search requires restarts" };

            unsigned n_threads = how_many_threads(params.n_threads);
            ThreadedSolver solver(model, params, n_threads);
            result = solver.solve();
        }

        if (params.proof && result.complete && result.mapping.empty())
            params.proof->finish_unsat_proof();

        return result;
    }
}

