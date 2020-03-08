/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "common_subgraph.hh"
#include "proof.hh"
#include "configuration.hh"

#include <algorithm>
#include <map>
#include <set>
#include <string_view>
#include <utility>
#include <vector>

using std::map;
using std::min;
using std::nullopt;
using std::pair;
using std::set;
using std::string;
using std::string_view;
using std::vector;

namespace
{
    enum class SearchResult
    {
        Aborted,
        Complete,
        DecidedTrue
    };

    struct SplitDomains
    {
        vector<pair<set<int>, set<int> > > partitions;
    };

    struct Assignments
    {
        vector<pair<int, int> > assigned;
        vector<int> rejected;
    };

    auto assignments_as_proof_decisions(const Assignments & assignments) -> vector<pair<int, int> >
    {
        vector<pair<int, int> > trail;
        for (auto & [ l, r ] : assignments.assigned)
            trail.emplace_back(l, r);
        return trail;
    }

    struct CommonSubgraphRunner
    {
        const InputGraph & first;
        const InputGraph & second;
        const CommonSubgraphParams & params;

        CommonSubgraphRunner(const InputGraph & f, const InputGraph & s, const CommonSubgraphParams & p) :
            first(f),
            second(s),
            params(p)
        {
        }

        auto branch_assigning(const SplitDomains & d, int left, int right) -> SplitDomains
        {
            SplitDomains result;

            for (auto & [ l, r ] : d.partitions) {
                set<int> new_l_1, new_r_1, new_l_2, new_r_2;

                for (auto & v : l) {
                    if (left != v) {
                        if (first.adjacent(left, v))
                            new_l_1.emplace(v);
                        else
                            new_l_2.emplace(v);
                    }
                }

                for (auto & v : r) {
                    if (right != v) {
                        if (second.adjacent(right, v))
                            new_r_1.emplace(v);
                        else
                            new_r_2.emplace(v);
                    }
                }

                if ((! new_l_1.empty()) && (! new_r_1.empty()))
                    result.partitions.emplace_back(new_l_1, new_r_1);
                if ((! new_l_2.empty()) && (! new_r_2.empty()))
                    result.partitions.emplace_back(new_l_2, new_r_2);
            }

            return result;
        }

        auto branch_rejecting(const SplitDomains & d, int left) -> SplitDomains
        {
            SplitDomains result;

            for (auto & [ l, r ] : d.partitions) {
                if (! l.count(left))
                    result.partitions.emplace_back(l, r);
                else {
                    set<int> new_l = l;
                    new_l.erase(left);
                    if (! new_l.empty())
                        result.partitions.emplace_back(new_l, r);
                }
            }

            return result;
        }

        auto bound(const SplitDomains & d) -> unsigned
        {
            unsigned result = 0;
            for (auto & [ l, r ] : d.partitions)
                result += min(l.size(), r.size());
            return result;
        };

        auto find_branch_partition(const SplitDomains & d) -> vector<pair<set<int>, set<int> > >::const_iterator
        {
            auto result = d.partitions.begin();

            for (auto b = d.partitions.begin(), b_end = d.partitions.end() ; b != b_end ; ++b)
                if (b->first.size() < result->first.size())
                    result = b;

            return result;
        };

        auto search(
                int depth,
                Assignments & assignments,
                Assignments & incumbent,
                const SplitDomains & domains,
                unsigned long long & nodes) -> SearchResult
        {
            ++nodes;

            auto branch_domains = find_branch_partition(domains);
            if (branch_domains == domains.partitions.end()) {
                if (assignments.assigned.size() > incumbent.assigned.size()) {
                    if (params.decide) {
                       if (assignments.assigned.size() >= *params.decide) {
                           incumbent = assignments;
                           return SearchResult::DecidedTrue;
                       }
                    }
                    else
                        incumbent = assignments;
                }
            }
            else {
                int left_branch = *branch_domains->first.begin();
                for (auto & right_branch : branch_domains->second) {
                    // branch with left_branch assigned to right_branch
                    if (params.proof) {
                        params.proof->guessing(depth, NamedVertex{ left_branch, first.vertex_name(left_branch) },
                                NamedVertex{ right_branch, second.vertex_name(right_branch) });
                        params.proof->start_level(depth + 1);
                    }

                    auto new_domains = branch_assigning(domains, left_branch, right_branch);
                    assignments.assigned.emplace_back(left_branch, right_branch);
                    if (assignments.assigned.size() + bound(new_domains) > incumbent.assigned.size()) {
                        switch (search(depth + 1, assignments, incumbent, new_domains, nodes)) {
                            case SearchResult::Aborted:     return SearchResult::Aborted;
                            case SearchResult::DecidedTrue: return SearchResult::DecidedTrue;
                            case SearchResult::Complete:    break;
                        }
                    }
                    else if (params.proof)
                        params.proof->mcs_bound(new_domains.partitions);

                    if (params.proof) {
                        params.proof->start_level(depth);
                        params.proof->incorrect_guess(assignments_as_proof_decisions(assignments), true);
                        params.proof->forget_level(depth + 1);
                    }

                    assignments.assigned.pop_back();
                }

                // now with left_branch assigned to null
                if (params.proof) {
                    params.proof->guessing(depth, NamedVertex{ left_branch, first.vertex_name(left_branch) },
                            NamedVertex{ second.size(), "null" });
                    params.proof->start_level(depth + 1);
                }

                auto new_domains = branch_rejecting(domains, left_branch);
                assignments.rejected.emplace_back(left_branch);
                if (assignments.assigned.size() + bound(new_domains) > incumbent.assigned.size()) {
                    switch (search(depth + 1, assignments, incumbent, new_domains, nodes)) {
                        case SearchResult::Aborted:     return SearchResult::Aborted;
                        case SearchResult::DecidedTrue: return SearchResult::DecidedTrue;
                        case SearchResult::Complete:    break;
                    }
                }
                else if (params.proof)
                    params.proof->mcs_bound(new_domains.partitions);
                if (params.proof) {
                    params.proof->start_level(depth);
                    params.proof->incorrect_guess(assignments_as_proof_decisions(assignments), true);
                    params.proof->forget_level(depth + 1);
                }
                assignments.rejected.pop_back();
            }

            if (params.proof)
                params.proof->out_of_guesses(assignments_as_proof_decisions(assignments));

            return SearchResult::Complete;
        }

        auto run() -> CommonSubgraphResult
        {
            CommonSubgraphResult result;

            map<pair<bool, string_view>, pair<set<int>, set<int> > > initial_partitions;

            for (int v = 0 ; v < first.size() ; ++v)
                initial_partitions[pair{ first.adjacent(v, v), first.vertex_label(v) }].first.insert(v);
            for (int v = 0 ; v < second.size() ; ++v)
                initial_partitions[pair{ second.adjacent(v, v), second.vertex_label(v) }].second.insert(v);

            SplitDomains domains;
            for (auto & [ k, p ] : initial_partitions) {
                auto & [ l, r ] = p;
                if ((! l.empty()) && (! r.empty()))
                    domains.partitions.emplace_back(l, r);
            }

            Assignments assignments, incumbent;

            if (params.decide)
                for (unsigned i = 1 ; i <= *params.decide - 1 ; ++i)
                    incumbent.assigned.emplace_back(-int(i), -int(i));

            if (params.decide && (bound(domains) < *params.decide)) {
                result.complete = true;
                if (params.proof)
                    params.proof->mcs_bound(domains.partitions);
            }
            else {
                switch (search(0, assignments, incumbent, domains, result.nodes)) {
                    case SearchResult::Aborted:
                        break;

                    case SearchResult::DecidedTrue:
                        result.complete = true;
                        for (auto & [ f, s ] : incumbent.assigned)
                            result.mapping.emplace(f, s);
                        break;

                    case SearchResult::Complete:
                        result.complete = true;
                        if (! params.decide) {
                            for (auto & [ f, s ] : incumbent.assigned)
                                result.mapping.emplace(f, s);
                        }
                        break;
                }
            }

            if (params.proof && result.complete && result.mapping.empty())
                params.proof->finish_unsat_proof();

            return result;
        }
    };
}

auto solve_common_subgraph_problem(const InputGraph & first, const InputGraph & second, const CommonSubgraphParams & params) -> CommonSubgraphResult
{
    if (params.proof) {
        if (! params.decide)
            throw UnsupportedConfiguration{ "Proof logging currently only works with decision problems" };

        for (int n = 0 ; n < first.size() ; ++n) {
            params.proof->create_cp_variable(n, second.size() + 1,
                    [&] (int v) { return first.vertex_name(v); },
                    [&] (int v) { if (v == second.size()) return string("null"); else return second.vertex_name(v); });
        }

        // generate constraints for injectivity
        params.proof->create_injectivity_constraints(first.size(), second.size());

        // generate edge constraints, and also handle loops here
        for (int p = 0 ; p < first.size() ; ++p) {
            for (int t = 0 ; t < second.size() ; ++t) {
                if (first.adjacent(p, p) && ! second.adjacent(t, t))
                    params.proof->create_forbidden_assignment_constraint(p, t);
                else if (first.vertex_label(p) != second.vertex_label(t))
                    params.proof->create_forbidden_assignment_constraint(p, t);
                else {
                    params.proof->start_adjacency_constraints_for(p, t);

                    // if p can be mapped to t, then each (non-)neighbour of p...
                    for (int q = 0 ; q < first.size() ; ++q)
                        if (q != p) {
                            // ... must be mapped to a (non-)neighbour of t
                            vector<int> permitted;
                            for (int u = 0 ; u < second.size() ; ++u)
                                if (t != u && first.adjacent(p, q) == second.adjacent(t, u))
                                    permitted.push_back(u);
                            // or null
                            permitted.push_back(second.size());
                            params.proof->create_adjacency_constraint(p, q, t, permitted);
                        }
                }
            }
        }

        params.proof->create_non_null_decision_bound(first.size(), second.size(), *params.decide);

        // output the model file
        params.proof->finalise_model();

        // rewrite the objective to the form we need
        params.proof->rewrite_mcs_objective(first.size());
    }

    CommonSubgraphRunner runner{ first, second, params };
    return runner.run();
}

