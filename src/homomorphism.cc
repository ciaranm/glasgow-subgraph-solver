/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "homomorphism.hh"
#include "cheap_all_different.hh"
#include "clique.hh"
#include "configuration.hh"
#include "graph_traits.hh"
#include "homomorphism_domain.hh"
#include "homomorphism_traits.hh"
#include "thread_utils.hh"
#include "watches.hh"
#include "proof.hh"
#include "svo_bitset.hh"
#include "subgraph_model.hh"

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <random>
#include <thread>
#include <tuple>
#include <type_traits>
#include <utility>

#include <boost/thread/barrier.hpp>

using std::atomic;
using std::conditional_t;
using std::fill;
using std::find_if;
using std::function;
using std::greater;
using std::is_same;
using std::make_optional;
using std::make_unique;
using std::max;
using std::map;
using std::move;
using std::mt19937;
using std::mutex;
using std::numeric_limits;
using std::optional;
using std::pair;
using std::sort;
using std::stable_sort;
using std::string;
using std::swap;
using std::thread;
using std::to_string;
using std::tuple;
using std::uniform_int_distribution;
using std::unique_lock;
using std::unique_ptr;
using std::vector;

using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::chrono::steady_clock;
using std::chrono::operator""ms;

using boost::barrier;

namespace
{
    enum class SearchResult
    {
        Aborted,
        Unsatisfiable,
        Satisfiable,
        SatisfiableButKeepGoing,
        Restart
    };

    struct Assignment
    {
        unsigned pattern_vertex;
        unsigned target_vertex;

        auto operator== (const Assignment & other) const -> bool
        {
            return pattern_vertex == other.pattern_vertex && target_vertex == other.target_vertex;
        }

        auto operator!= (const Assignment & other) const -> bool
        {
            return ! (*this == other);
        }
    };

    struct AssignmentInformation
    {
        Assignment assignment;
        bool is_decision;
        int discrepancy_count;
        int choice_count;
    };

    struct Assignments
    {
        vector<AssignmentInformation> values;

        bool contains(const Assignment & assignment) const
        {
            // this should not be a linear scan...
            return values.end() != find_if(values.begin(), values.end(), [&] (const auto & a) {
                    return a.assignment == assignment;
                    });
        }
    };

    template <typename EntryType_>
    struct AssignmentWatchTable
    {
        unsigned target_size;
        vector<EntryType_> data;

        EntryType_ & operator[] (Assignment x)
        {
            return data[target_size * x.pattern_vertex + x.target_vertex];
        }
    };

    struct Searcher
    {
        using Domains = vector<HomomorphismDomain>;

        const SubgraphModel & model;
        const HomomorphismParams & params;

        Watches<Assignment, AssignmentWatchTable> watches;

        mt19937 global_rand;

        Searcher(const SubgraphModel & m, const HomomorphismParams & p) :
            model(m),
            params(p)
        {
            // set up space for watches
            if (might_have_watches(params)) {
                watches.table.target_size = model.target_size;
                watches.table.data.resize(model.pattern_size * model.target_size);
            }
        }

        auto assignments_as_proof_decisions(const Assignments & assignments) const -> vector<pair<int, int> >
        {
            vector<pair<int, int> > trail;
            for (auto & a : assignments.values)
                if (a.is_decision)
                    trail.emplace_back(a.assignment.pattern_vertex, a.assignment.target_vertex);
            return trail;
        }

        auto solution_in_proof_form(const Assignments & assignments) const -> vector<pair<NamedVertex, NamedVertex> >
        {
            vector<pair<NamedVertex, NamedVertex> > solution;
            for (auto & a : assignments.values)
                if (solution.end() == find_if(solution.begin(), solution.end(),
                            [&] (const auto & t) { return unsigned(t.first.first) == a.assignment.pattern_vertex; }))
                    solution.emplace_back(
                            model.pattern_vertex_for_proof(a.assignment.pattern_vertex),
                            model.target_vertex_for_proof(a.assignment.target_vertex));
            return solution;
        }

        auto expand_to_full_result(const Assignments & assignments, VertexToVertexMapping & mapping) -> void
        {
            for (auto & a : assignments.values)
                mapping.emplace(a.assignment.pattern_vertex, a.assignment.target_vertex);
        }

        auto find_unit_domain(Domains & domains) -> typename Domains::iterator
        {
            return find_if(domains.begin(), domains.end(), [] (HomomorphismDomain & d) {
                    return (! d.fixed) && 1 == d.count;
                    });
        }

        template <bool has_edge_labels_, bool induced_>
        auto propagate_adjacency_constraints(HomomorphismDomain & d, const Assignment & current_assignment) -> void
        {
            auto graph_pairs_to_consider = model.pattern_adjacency_bits(current_assignment.pattern_vertex, d.v);

            // for the original graph pair, if we're adjacent...
            if (graph_pairs_to_consider & (1u << 0)) {
                // ...then we can only be mapped to adjacent vertices
                d.values &= model.target_graph_row(0, current_assignment.target_vertex);
            }
            else {
                if constexpr (induced_) {
                    // ...otherwise we can only be mapped to adjacent vertices
                    d.values.intersect_with_complement(model.target_graph_row(0, current_assignment.target_vertex));
                }
            }

            // and for each remaining graph pair...
            for (unsigned g = 1 ; g < model.max_graphs ; ++g) {
                // if we're adjacent...
                if (graph_pairs_to_consider & (1u << g)) {
                    // ...then we can only be mapped to adjacent vertices
                    d.values &= model.target_graph_row(g, current_assignment.target_vertex);
                }
            }

            if constexpr (has_edge_labels_) {
                // if we're adjacent in the original graph, additionally the edge labels need to match up
                if (graph_pairs_to_consider & (1u << 0)) {
                    auto check_d_values = d.values;

                    auto want_forward_label = model.pattern_edge_label(d.v, current_assignment.pattern_vertex);
                    auto want_reverse_label = model.pattern_edge_label(current_assignment.pattern_vertex, d.v);
                    for (auto c = check_d_values.find_first() ; c != decltype(check_d_values)::npos ; c = check_d_values.find_first()) {
                        check_d_values.reset(c);

                        auto got_forward_label = model.target_edge_label(c, current_assignment.target_vertex);
                        auto got_reverse_label = model.target_edge_label(current_assignment.target_vertex, c);

                        if (got_forward_label != want_forward_label || got_reverse_label != want_reverse_label)
                            d.values.reset(c);
                    }
                }
            }
        }

        auto both_in_the_neighbourhood_of_some_vertex(unsigned v, unsigned w) -> bool
        {
            auto i = model.pattern_graph_row(0, v);
            i &= model.pattern_graph_row(0, w);
            return i.any();
        }

        auto propagate_simple_constraints(Domains & new_domains, const Assignment & current_assignment) -> bool
        {
            // propagate for each remaining domain...
            for (auto & d : new_domains) {
                if (d.fixed)
                    continue;

                // injectivity
                switch (params.injectivity) {
                    case Injectivity::Injective:
                        d.values.reset(current_assignment.target_vertex);
                        break;
                    case Injectivity::LocallyInjective:
                        if (both_in_the_neighbourhood_of_some_vertex(current_assignment.pattern_vertex, d.v))
                            d.values.reset(current_assignment.target_vertex);
                        break;
                    case Injectivity::NonInjective:
                        break;
                }

                // adjacency
                if (! model.has_edge_labels()) {
                    if (params.induced)
                        propagate_adjacency_constraints<false, true>(d, current_assignment);
                    else
                        propagate_adjacency_constraints<false, false>(d, current_assignment);
                }
                else {
                    if (params.induced)
                        propagate_adjacency_constraints<true, true>(d, current_assignment);
                    else
                        propagate_adjacency_constraints<true, false>(d, current_assignment);
                    break;
                }

                // we might have removed values
                d.count = d.values.count();
                if (0 == d.count)
                    return false;
            }

            return true;
        }

        auto propagate_less_thans(Domains & new_domains) -> bool
        {
            vector<int> find_domain(model.pattern_size, -1);

            for (unsigned i = 0, i_end = new_domains.size() ; i != i_end ; ++i)
                find_domain[new_domains[i].v] = i;

            for (auto & [ a, b ] : model.pattern_less_thans_in_convenient_order) {
                if (find_domain[a] == -1 || find_domain[b] == -1)
                    continue;
                auto & a_domain = new_domains[find_domain[a]];
                auto & b_domain = new_domains[find_domain[b]];

               // first value of b must be at least one after the first possible value of a
               auto first_a = a_domain.values.find_first();
               if (first_a == decltype(a_domain.values)::npos)
                   return false;
               auto first_allowed_b = first_a + 1;

               if (first_allowed_b >= model.target_size)
                   return false;

               for (auto v = b_domain.values.find_first() ; v != decltype(b_domain.values)::npos ; v = b_domain.values.find_first()) {
                   if (v >= first_allowed_b)
                       break;
                   b_domain.values.reset(v);
               }

               // b might have shrunk (and detect empty before the next bit to make life easier)
               b_domain.count = b_domain.values.count();
               if (0 == b_domain.count)
                   return false;
            }

            for (auto & [ a, b ] : model.pattern_less_thans_in_convenient_order) {
                if (find_domain[a] == -1 || find_domain[b] == -1)
                    continue;
                auto & a_domain = new_domains[find_domain[a]];
                auto & b_domain = new_domains[find_domain[b]];

                // last value of a must be at least one before the last possible value of b
                auto b_values_copy = b_domain.values;
                auto last_b = b_domain.values.find_first();
                for (auto v = last_b ; v != decltype(b_values_copy)::npos ; v = b_values_copy.find_first()) {
                    b_values_copy.reset(v);
                    last_b = v;
                }

                if (last_b == 0)
                    return false;
                auto last_allowed_a = last_b - 1;

                auto a_values_copy = a_domain.values;
                for (auto v = a_values_copy.find_first() ; v != decltype(a_values_copy)::npos ; v = a_values_copy.find_first()) {
                    a_values_copy.reset(v);
                    if (v > last_allowed_a)
                        a_domain.values.reset(v);
                }

                // a might have shrunk
                a_domain.count = a_domain.values.count();
                if (0 == a_domain.count)
                    return false;
            }

            return true;
        }

        auto propagate(Domains & new_domains, Assignments & assignments) -> bool
        {
            // whilst we've got a unit domain...
            for (typename Domains::iterator branch_domain = find_unit_domain(new_domains) ;
                    branch_domain != new_domains.end() ;
                    branch_domain = find_unit_domain(new_domains)) {
                // what are we assigning?
                Assignment current_assignment = { branch_domain->v, unsigned(branch_domain->values.find_first()) };

                // ok, make the assignment
                branch_domain->fixed = true;
                assignments.values.push_back({ current_assignment, false, -1, -1 });

                if (params.proof)
                    params.proof->unit_propagating(
                            model.pattern_vertex_for_proof(current_assignment.pattern_vertex),
                            model.target_vertex_for_proof(current_assignment.target_vertex));

                // propagate watches
                if (might_have_watches(params))
                    watches.propagate(current_assignment,
                            [&] (const Assignment & a) { return ! assignments.contains(a); },
                            [&] (const Assignment & a) {
                                    for (auto & d : new_domains) {
                                        if (d.fixed)
                                            continue;

                                        if (d.v == a.pattern_vertex) {
                                            d.values.reset(a.target_vertex);
                                            break;
                                        }
                                    }
                                });

                // propagate simple all different and adjacency
                if (! propagate_simple_constraints(new_domains, current_assignment))
                    return false;

                // propagate less than
                if (model.has_less_thans() && ! propagate_less_thans(new_domains))
                    return false;

                // propagate all different
                if (params.injectivity == Injectivity::Injective)
                    if (! cheap_all_different(model.target_size, new_domains, params.proof))
                        return false;
            }

            return true;
        }

        auto find_branch_domain(const Domains & domains) -> const HomomorphismDomain *
        {
            const HomomorphismDomain * result = nullptr;
            for (auto & d : domains)
                if (! d.fixed)
                    if ((! result) ||
                            (d.count < result->count) ||
                            (d.count == result->count && model.pattern_degree(0, d.v) > model.pattern_degree(0, result->v)))
                        result = &d;
            return result;
        }

        auto copy_nonfixed_domains_and_make_assignment(
                const Domains & domains,
                unsigned branch_v,
                unsigned f_v) -> Domains
        {
            Domains new_domains;
            new_domains.reserve(domains.size());
            for (auto & d : domains) {
                if (d.fixed)
                    continue;

                new_domains.push_back(d);
                if (d.v == branch_v) {
                    new_domains.back().values.reset();
                    new_domains.back().values.set(f_v);
                    new_domains.back().count = 1;
                }
            }
            return new_domains;
        }

        auto post_nogood(
                const Assignments & assignments)
        {
            if (! might_have_watches(params))
                return;

            Nogood<Assignment> nogood;

            for (auto & a : assignments.values)
                if (a.is_decision)
                    nogood.literals.emplace_back(a.assignment);

            watches.post_nogood(move(nogood));

            if (params.proof)
                params.proof->post_restart_nogood(assignments_as_proof_decisions(assignments));
        }

        auto softmax_shuffle(
                vector<int> & branch_v,
                unsigned branch_v_end
                ) -> void
        {
            // repeatedly pick a softmax-biased vertex, move it to the front of branch_v,
            // and then only consider items further to the right in the next iteration.

            // Using floating point here turned out to be way too slow. Fortunately the base
            // of softmax doesn't seem to matter, so we use 2 instead of e, and do everything
            // using bit voodoo.
            auto expish = [largest_target_degree = this->model.largest_target_degree()] (int degree) {
                constexpr int sufficient_space_for_adding_up = numeric_limits<long long>::digits - 18;
                auto shift = max<int>(degree - largest_target_degree + sufficient_space_for_adding_up, 0);
                return 1ll << shift;
            };

            long long total = 0;
            for (unsigned v = 0 ; v < branch_v_end ; ++v)
                total += expish(model.target_degree(0, branch_v[v]));

            for (unsigned start = 0 ; start < branch_v_end ; ++start) {
                // pick a random number between 1 and total inclusive
                uniform_int_distribution<long long> dist(1, total);
                long long select_score = dist(global_rand);

                // go over the list until we hit the score
                unsigned select_element = start;
                for ( ; select_element + 1 < branch_v_end ; ++select_element) {
                    select_score -= expish(model.target_degree(0, branch_v[select_element]));
                    if (select_score <= 0)
                        break;
                }

                // move to front
                total -= expish(model.target_degree(0, branch_v[select_element]));
                swap(branch_v[select_element], branch_v[start]);
            }
        }

        auto degree_sort(
                vector<int> & branch_v,
                unsigned branch_v_end,
                bool reverse
                ) -> void
        {
            stable_sort(branch_v.begin(), branch_v.begin() + branch_v_end, [&] (int a, int b) -> bool {
                    return (model.target_degree(0, a) >= model.target_degree(0, b)) ^ reverse;
                    });
        }

        auto restarting_search(
                Assignments & assignments,
                const Domains & domains,
                unsigned long long & nodes,
                unsigned long long & propagations,
                unsigned long long & solution_count,
                int depth,
                RestartsSchedule & restarts_schedule) -> SearchResult
        {
            if (params.timeout->should_abort())
                return SearchResult::Aborted;

            ++nodes;

            // find ourselves a domain, or succeed if we're all assigned
            const HomomorphismDomain * branch_domain = find_branch_domain(domains);
            if (! branch_domain) {
                if (params.lackey) {
                    VertexToVertexMapping mapping;
                    expand_to_full_result(assignments, mapping);
                    if (! params.lackey->check_solution(mapping))
                        return SearchResult::Unsatisfiable;
                }

                if (params.proof)
                    params.proof->post_solution(solution_in_proof_form(assignments));

                if (params.count_solutions) {
                    ++solution_count;
                    if (params.enumerate_callback) {
                        VertexToVertexMapping mapping;
                        expand_to_full_result(assignments, mapping);
                        params.enumerate_callback(mapping);
                    }

                    return SearchResult::SatisfiableButKeepGoing;
                }
                else
                    return SearchResult::Satisfiable;
            }

            // pull out the remaining values in this domain for branching
            auto remaining = branch_domain->values;

            vector<int> branch_v(model.target_size);

            unsigned branch_v_end = 0;
            for (auto f_v = remaining.find_first() ; f_v != decltype(remaining)::npos ; f_v = remaining.find_first()) {
                remaining.reset(f_v);
                branch_v[branch_v_end++] = f_v;
            }

            switch (params.value_ordering_heuristic) {
                case ValueOrdering::Degree:
                    degree_sort(branch_v, branch_v_end, false);
                    break;

                case ValueOrdering::AntiDegree:
                    degree_sort(branch_v, branch_v_end, true);
                    break;

                case ValueOrdering::Biased:
                    softmax_shuffle(branch_v, branch_v_end);
                    break;

                case ValueOrdering::Random:
                    shuffle(branch_v.begin(), branch_v.begin() + branch_v_end, global_rand);
                    break;
            }

            int discrepancy_count = 0;
            bool actually_hit_a_failure = false;

            // for each value remaining...
            for (auto f_v = branch_v.begin(), f_end = branch_v.begin() + branch_v_end ; f_v != f_end ; ++f_v) {
                if (params.proof)
                    params.proof->guessing(depth, model.pattern_vertex_for_proof(branch_domain->v), model.target_vertex_for_proof(*f_v));

                // modified in-place by appending, we can restore by shrinking
                auto assignments_size = assignments.values.size();

                // make the assignment
                assignments.values.push_back({ { branch_domain->v, unsigned(*f_v) }, true, discrepancy_count, int(branch_v_end) });

                // set up new domains
                Domains new_domains = copy_nonfixed_domains_and_make_assignment(domains, branch_domain->v, *f_v);

                // propagate
                ++propagations;
                if (! propagate(new_domains, assignments)) {
                    // failure? restore assignments and go on to the next thing
                    if (params.proof)
                        params.proof->propagation_failure(assignments_as_proof_decisions(assignments), model.pattern_vertex_for_proof(branch_domain->v), model.target_vertex_for_proof(*f_v));

                    assignments.values.resize(assignments_size);
                    actually_hit_a_failure = true;

                    continue;
                }

                if (params.proof)
                    params.proof->start_level(depth + 2);

                // recursive search
                auto search_result = restarting_search(assignments, new_domains, nodes, propagations,
                        solution_count, depth + 1, restarts_schedule);

                switch (search_result) {
                    case SearchResult::Satisfiable:
                        return SearchResult::Satisfiable;

                    case SearchResult::Aborted:
                        return SearchResult::Aborted;

                    case SearchResult::Restart:
                        // restore assignments before posting nogoods, it's easier
                        assignments.values.resize(assignments_size);

                        // post nogoods for everything we've done so far
                        for (auto l = branch_v.begin() ; l != f_v ; ++l) {
                            assignments.values.push_back({ { branch_domain->v, unsigned(*l) }, true, -2, -2 });
                            post_nogood(assignments);
                            assignments.values.pop_back();
                        }

                        return SearchResult::Restart;

                    case SearchResult::SatisfiableButKeepGoing:
                        if (params.proof) {
                            params.proof->back_up_to_level(depth + 1);
                            params.proof->incorrect_guess(assignments_as_proof_decisions(assignments), false);
                            params.proof->forget_level(depth + 2);
                        }

                        // restore assignments
                        assignments.values.resize(assignments_size);
                        break;

                    case SearchResult::Unsatisfiable:
                        if (params.proof) {
                            params.proof->back_up_to_level(depth + 1);
                            params.proof->incorrect_guess(assignments_as_proof_decisions(assignments), true);
                            params.proof->forget_level(depth + 2);
                        }

                        // restore assignments
                        assignments.values.resize(assignments_size);
                        actually_hit_a_failure = true;

                        break;
                }

                ++discrepancy_count;
            }

            // no values remaining, backtrack, or possibly kick off a restart
            if (params.proof)
                params.proof->out_of_guesses(assignments_as_proof_decisions(assignments));

            if (actually_hit_a_failure)
                restarts_schedule.did_a_backtrack();

            if (restarts_schedule.should_restart()) {
                if (params.proof)
                    params.proof->back_up_to_top();
                post_nogood(assignments);
                return SearchResult::Restart;
            }
            else
                return SearchResult::Unsatisfiable;
        }

        auto save_result(const Assignments & assignments, HomomorphismResult & result) -> void
        {
            expand_to_full_result(assignments, result.mapping);

            string where = "where =";
            for (auto & a : assignments.values)
                where.append(" " + to_string(a.discrepancy_count) + "/" + to_string(a.choice_count));
            result.extra_stats.push_back(where);
        }
    };

    struct HomomorphismSolver
    {
        using Domains = vector<HomomorphismDomain>;

        const SubgraphModel & model;
        const HomomorphismParams & params;

        HomomorphismSolver(const SubgraphModel & m, const HomomorphismParams & p) :
            model(m),
            params(p)
        {
        }

        auto check_label_compatibility(int p, int t) -> bool
        {
            if (! model.has_vertex_labels())
                return true;
            else
                return model.pattern_vertex_label(p) == model.target_vertex_label(t);
        }

        auto check_loop_compatibility(int p, int t) -> bool
        {
            if (model.pattern_has_loop(p) && ! model.target_has_loop(t))
                return false;
            else if (params.induced && (model.pattern_has_loop(p) != model.target_has_loop(t)))
                return false;

            return true;
        }

        auto check_degree_compatibility(
                int p,
                int t,
                unsigned graphs_to_consider,
                vector<vector<vector<int> > > & patterns_ndss,
                vector<vector<optional<vector<int> > > > & targets_ndss,
                bool do_not_do_nds_yet
                ) -> bool
        {
            if (! degree_and_nds_are_preserved(params))
                return true;

            for (unsigned g = 0 ; g < graphs_to_consider ; ++g) {
                if (model.target_degree(g, t) < model.pattern_degree(g, p)) {
                    // not ok, degrees differ
                    if (params.proof) {
                        // get the actual neighbours of p and t, in their original terms
                        vector<int> n_p, n_t;

                        auto np = model.pattern_graph_row(g, p);
                        for (unsigned j = 0 ; j < model.pattern_size ; ++j)
                            if (np.test(j))
                                n_p.push_back(j);

                        auto nt = model.target_graph_row(g, t);
                        for (auto j = nt.find_first() ; j != decltype(nt)::npos ; j = nt.find_first()) {
                            nt.reset(j);
                            n_t.push_back(j);
                        }

                        params.proof->incompatible_by_degrees(g, model.pattern_vertex_for_proof(p), n_p,
                                model.target_vertex_for_proof(t), n_t);
                    }
                    return false;
                }
                else if (degree_and_nds_are_exact(params, model.pattern_size, model.target_size)
                        && model.target_degree(g, t) != model.pattern_degree(g, p)) {
                    // not ok, degrees must be exactly the same
                    return false;
                }
            }
            if (params.no_nds || do_not_do_nds_yet)
                return true;

            // full compare of neighbourhood degree sequences
            if (! targets_ndss.at(0).at(t)) {
                for (unsigned g = 0 ; g < graphs_to_consider ; ++g) {
                    targets_ndss.at(g).at(t) = vector<int>{};
                    auto ni = model.target_graph_row(g, t);
                    for (auto j = ni.find_first() ; j != decltype(ni)::npos ; j = ni.find_first()) {
                        ni.reset(j);
                        targets_ndss.at(g).at(t)->push_back(model.target_degree(g, j));
                    }
                    sort(targets_ndss.at(g).at(t)->begin(), targets_ndss.at(g).at(t)->end(), greater<int>());
                }
            }

            for (unsigned g = 0 ; g < graphs_to_consider ; ++g) {
                for (unsigned x = 0 ; x < patterns_ndss.at(g).at(p).size() ; ++x) {
                    if (targets_ndss.at(g).at(t)->at(x) < patterns_ndss.at(g).at(p).at(x)) {
                        if (params.proof) {
                            vector<int> p_subsequence, t_subsequence, t_remaining;

                            // need to know the NDS together with the actual vertices
                            vector<pair<int, int> > p_nds, t_nds;

                            auto np = model.pattern_graph_row(g, p);
                            for (auto w = np.find_first() ; w != decltype(np)::npos ; w = np.find_first()) {
                                np.reset(w);
                                p_nds.emplace_back(w, model.pattern_graph_row(g, w).count());
                            }

                            auto nt = model.target_graph_row(g, t);
                            for (auto w = nt.find_first() ; w != decltype(nt)::npos ; w = nt.find_first()) {
                                nt.reset(w);
                                t_nds.emplace_back(w, model.target_graph_row(g, w).count());
                            }

                            sort(p_nds.begin(), p_nds.end(), [] (const pair<int, int> & a, const pair<int, int> & b) {
                                    return a.second > b.second; });
                            sort(t_nds.begin(), t_nds.end(), [] (const pair<int, int> & a, const pair<int, int> & b) {
                                    return a.second > b.second; });

                            for (unsigned y = 0 ; y <= x ; ++y) {
                                p_subsequence.push_back(p_nds[y].first);
                                t_subsequence.push_back(t_nds[y].first);
                            }
                            for (unsigned y = x + 1 ; y < t_nds.size() ; ++y)
                                t_remaining.push_back(t_nds[y].first);

                            params.proof->incompatible_by_nds(g, model.pattern_vertex_for_proof(p),
                                    model.target_vertex_for_proof(t), p_subsequence, t_subsequence, t_remaining);
                        }
                        return false;
                    }
                    else if (degree_and_nds_are_exact(params, model.pattern_size, model.target_size)
                            && targets_ndss.at(g).at(t)->at(x) != patterns_ndss.at(g).at(p).at(x))
                        return false;
                }
            }

            return true;
        }

        auto initialise_domains(Domains & domains) -> bool
        {
            unsigned graphs_to_consider = model.max_graphs;

            /* pattern and target neighbourhood degree sequences */
            vector<vector<vector<int> > > patterns_ndss(graphs_to_consider);
            vector<vector<optional<vector<int> > > > targets_ndss(graphs_to_consider);

            if (degree_and_nds_are_preserved(params) && ! params.no_nds) {
                for (unsigned g = 0 ; g < graphs_to_consider ; ++g) {
                    patterns_ndss.at(g).resize(model.pattern_size);
                    targets_ndss.at(g).resize(model.target_size);
                }

                for (unsigned g = 0 ; g < graphs_to_consider ; ++g) {
                    for (unsigned i = 0 ; i < model.pattern_size ; ++i) {
                        auto ni = model.pattern_graph_row(g, i);
                        for (auto j = ni.find_first() ; j != decltype(ni)::npos ; j = ni.find_first()) {
                            ni.reset(j);
                            patterns_ndss.at(g).at(i).push_back(model.pattern_degree(g, j));
                        }
                        sort(patterns_ndss.at(g).at(i).begin(), patterns_ndss.at(g).at(i).end(), greater<int>());
                    }
                }
            }

            for (unsigned i = 0 ; i < model.pattern_size ; ++i) {
                domains.at(i).v = i;
                domains.at(i).values.reset();

                for (unsigned j = 0 ; j < model.target_size ; ++j) {
                    bool ok = true;

                    if (! check_label_compatibility(i, j))
                        ok = false;
                    else if (! check_loop_compatibility(i, j))
                        ok = false;
                    else if (! check_degree_compatibility(i, j, graphs_to_consider, patterns_ndss, targets_ndss, params.proof.get()))
                        ok = false;

                    if (ok)
                        domains.at(i).values.set(j);
                }

                domains.at(i).count = domains.at(i).values.count();
                if (0 == domains.at(i).count)
                    return false;
            }

            // for proof logging, we need degree information before we can output nds proofs
            if (params.proof && degree_and_nds_are_preserved(params) && ! params.no_nds) {
                for (unsigned i = 0 ; i < model.pattern_size ; ++i) {
                    for (unsigned j = 0 ; j < model.target_size ; ++j) {
                        if (domains.at(i).values.test(j) &&
                                ! check_degree_compatibility(i, j, graphs_to_consider, patterns_ndss, targets_ndss, false)) {
                            domains.at(i).values.reset(j);
                            if (0 == --domains.at(i).count)
                                return false;
                        }
                    }
                }
            }

            // quick sanity check that we have enough values
            if (is_nonshrinking(params)) {
                SVOBitset domains_union{ model.target_size, 0 };
                for (auto & d : domains)
                    domains_union |= d.values;

                unsigned domains_union_popcount = domains_union.count();
                if (domains_union_popcount < unsigned(model.pattern_size)) {
                    if (params.proof) {
                        vector<int> hall_lhs, hall_rhs;
                        for (auto & d : domains)
                            hall_lhs.push_back(d.v);
                        auto dd = domains_union;
                        for (auto v = dd.find_first() ; v != decltype(dd)::npos ; v = dd.find_first()) {
                            dd.reset(v);
                            hall_rhs.push_back(v);
                        }
                        params.proof->emit_hall_set_or_violator(hall_lhs, hall_rhs);
                    }
                    return false;
                }
            }

            for (auto & d : domains) {
                d.count = d.values.count();
                if (0 == d.count && params.proof) {
                    params.proof->initial_domain_is_empty(d.v);
                    return false;
                }
            }

            return true;
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
            if (! initialise_domains(domains)) {
                result.complete = true;
                return result;
            }

            // assignments
            Assignments assignments;
            assignments.values.reserve(model.pattern_size);

            // start search timer
            auto search_start_time = steady_clock::now();

            // do the search
            bool done = false;
            unsigned number_of_restarts = 0;

            Searcher searcher(model, params);

            while (! done) {
                ++number_of_restarts;

                // start watching new nogoods
                done = searcher.watches.apply_new_nogoods(
                        [&] (const Assignment & assignment) {
                            for (auto & d : domains)
                                if (d.v == assignment.pattern_vertex) {
                                    d.values.reset(assignment.target_vertex);
                                    d.count = d.values.count();
                                    break;
                                }
                        });

                if (done)
                    break;

                searcher.watches.clear_new_nogoods();

                ++result.propagations;
                if (searcher.propagate(domains, assignments)) {
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

            return result;
        }
    };

    struct ThreadedSolver : HomomorphismSolver
    {
        unsigned n_threads;

        ThreadedSolver(const SubgraphModel & m, const HomomorphismParams & p, unsigned t) :
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
            if (! initialise_domains(common_domains)) {
                common_result.complete = true;
                return common_result;
            }

            // start search timer
            auto search_start_time = steady_clock::now();

            vector<thread> threads;
            threads.reserve(n_threads);

            vector<unique_ptr<Searcher> > searchers{ n_threads };

            barrier wait_for_new_nogoods_barrier{ n_threads }, synced_nogoods_barrier{ n_threads };
            atomic<bool> restart_synchroniser{ false };

            function<auto (unsigned) -> void> work_function = [&searchers, &common_domains, &threads, &work_function,
                        &model = this->model, &params = this->params, n_threads = this->n_threads,
                        &common_result, &common_result_mutex, &by_thread_nodes, &by_thread_propagations,
                        &wait_for_new_nogoods_barrier, &synced_nogoods_barrier, &restart_synchroniser] (unsigned t) -> void
            {
                // do the search
                HomomorphismResult thread_result;

                bool just_the_first_thread = (0 == t) && params.delay_thread_creation;

                searchers[t] = make_unique<Searcher>(model, params);
                if (0 != t)
                    searchers[t]->global_rand.seed(t);

                unsigned number_of_restarts = 0;

                Domains domains = common_domains;

                Assignments thread_assignments;
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
                                [&] (const Assignment & assignment) {
                                    for (auto & d : domains)
                                        if (d.v == assignment.pattern_vertex) {
                                            d.values.reset(assignment.target_vertex);
                                            d.count = d.values.count();
                                            break;
                                        }
                                }))
                            break;

                        if (0 == t)
                            restart_synchroniser.store(false);

                        synced_nogoods_barrier.wait();

                        searchers[t]->watches.clear_new_nogoods();
                    }

                    ++thread_result.propagations;
                    if (searchers[t]->propagate(domains, thread_assignments)) {
                        auto assignments_copy = thread_assignments;

                        switch (searchers[t]->restarting_search(assignments_copy, domains, thread_result.nodes, thread_result.propagations,
                                    thread_result.solution_count, 0, *thread_restarts_schedule)) {
                            case SearchResult::Satisfiable:
                                searchers[t]->save_result(assignments_copy, thread_result);
                                thread_result.complete = true;
                                params.timeout->trigger_early_abort();
                                searchers[t]->watches.post_nogood(Nogood<Assignment>{ });
                                break;

                            case SearchResult::SatisfiableButKeepGoing:
                                thread_result.complete = true;
                                params.timeout->trigger_early_abort();
                                searchers[t]->watches.post_nogood(Nogood<Assignment>{ });
                                break;

                            case SearchResult::Unsatisfiable:
                                thread_result.complete = true;
                                params.timeout->trigger_early_abort();
                                searchers[t]->watches.post_nogood(Nogood<Assignment>{ });
                                break;

                            case SearchResult::Aborted:
                                searchers[t]->watches.post_nogood(Nogood<Assignment>{ });
                                break;

                            case SearchResult::Restart:
                                break;
                        }
                    }
                    else {
                        thread_result.complete = true;
                        searchers[t]->watches.post_nogood(Nogood<Assignment>{ });
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
        if (! params.pattern_less_constraints.empty())
            throw UnsupportedConfiguration{ "Proof logging cannot yet be used with less-constraints" };
        if (params.injectivity != Injectivity::Injective)
            throw UnsupportedConfiguration{ "Proof logging can currently only be used with injectivity" };
        if (params.induced)
            throw UnsupportedConfiguration{ "Proof logging cannot yet be used for induced problems" };
        if (pattern.has_vertex_labels() || pattern.has_edge_labels())
            throw UnsupportedConfiguration{ "Proof logging cannot yet be used on labelled graphs" };

        // set up our model file, with a set of OPB variables for each CP variable
        for (int n = 0 ; n < pattern.size() ; ++n) {
            params.proof->create_cp_variable(n, target.size(),
                    [&] (int v) { return pattern.vertex_name(v); },
                    [&] (int v) { return target.vertex_name(v); });
        }

        // generate constraints for injectivity
        params.proof->create_injectivity_constraints(pattern.size(), target.size());

        // generate edge constraints, and also handle loops here
        for (int p = 0 ; p < pattern.size() ; ++p) {
            for (int t = 0 ; t < target.size() ; ++t) {
                if (pattern.adjacent(p, p) && ! target.adjacent(t, t))
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
                        params.proof->create_adjacency_constraint(p, q, t, permitted);
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
        SubgraphModel model(target, pattern, params);

        if (! model.prepare(params)) {
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

