#include <gss/clique.hh>
#include <gss/common_subgraph.hh>
#include <gss/configuration.hh>
#include <gss/innards/proof.hh>

#include <algorithm>
#include <map>
#include <optional>
#include <set>
#include <string_view>
#include <utility>
#include <vector>

using namespace gss;
using namespace gss::innards;

using std::function;
using std::make_optional;
using std::make_shared;
using std::make_unique;
using std::map;
using std::min;
using std::nullopt;
using std::optional;
using std::pair;
using std::set;
using std::shared_ptr;
using std::string;
using std::string_view;
using std::tuple;
using std::vector;

namespace
{
    enum class SearchResult
    {
        Aborted,
        Complete,
        DecidedTrue,
        SatisfiableButKeepGoing
    };

    struct SplitDomains
    {
        vector<pair<set<int>, set<int>>> partitions;
    };

    struct Assignments
    {
        vector<pair<int, int>> assigned;
        vector<int> rejected;
    };

    auto assignments_as_proof_decisions(const Assignments & assignments) -> vector<pair<int, int>>
    {
        vector<pair<int, int>> trail;
        for (auto & [l, r] : assignments.assigned)
            trail.emplace_back(l, r);
        return trail;
    }

    struct CommonSubgraphRunner
    {
        const InputGraph & first;
        const InputGraph & second;
        const CommonSubgraphParams & params;
        shared_ptr<Proof> proof;

        CommonSubgraphRunner(const InputGraph & f, const InputGraph & s, const CommonSubgraphParams & p,
            const shared_ptr<Proof> & r) :
            first(f),
            second(s),
            params(p),
            proof(r)
        {
        }

        auto branch_assigning(const SplitDomains & d, int left, int right) -> SplitDomains
        {
            SplitDomains result;

            for (auto & [l, r] : d.partitions) {
                map<tuple<bool, bool, string_view, string_view>, pair<set<int>, set<int>>> new_partitions;

                string no_label;
                auto partition_of = [&](const InputGraph & g, int w, int v) -> tuple<bool, bool, string_view, string_view> {
                    return tuple{
                        g.adjacent(w, v),
                        g.adjacent(v, w),
                        g.adjacent(w, v) ? g.edge_label(w, v) : no_label,
                        g.adjacent(v, w) ? g.edge_label(v, w) : no_label};
                };

                for (auto & v : l)
                    if (left != v)
                        new_partitions[partition_of(first, left, v)].first.emplace(v);

                for (auto & v : r)
                    if (right != v)
                        new_partitions[partition_of(second, right, v)].second.emplace(v);

                for (auto & [_, s] : new_partitions)
                    if ((! s.first.empty()) && (! s.second.empty()))
                        result.partitions.emplace_back(s.first, s.second);
            }

            return result;
        }

        auto branch_rejecting(const SplitDomains & d, int left) -> SplitDomains
        {
            SplitDomains result;

            for (auto & [l, r] : d.partitions) {
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
            for (auto & [l, r] : d.partitions)
                result += min(l.size(), r.size());
            return result;
        };

        auto find_branch_partition(
            const SplitDomains & d,
            const optional<set<int>> & permitted_branch_variables) -> vector<pair<set<int>, set<int>>>::const_iterator
        {
            auto result = d.partitions.end();

            for (auto b = d.partitions.begin(), b_end = d.partitions.end(); b != b_end; ++b)
                if ((! permitted_branch_variables) || (permitted_branch_variables->count(*b->first.begin())))
                    if (result == d.partitions.end() || b->first.size() < result->first.size())
                        result = b;

            return result;
        };

        auto search(
            int depth,
            Assignments & assignments,
            Assignments & incumbent,
            const SplitDomains & domains,
            unsigned long long & nodes,
            loooong & solution_count,
            const optional<set<int>> & permitted_branch_variables) -> SearchResult
        {
            if (params.timeout->should_abort())
                return SearchResult::Aborted;

            ++nodes;

            auto branch_domains = find_branch_partition(domains, permitted_branch_variables);
            if (branch_domains == domains.partitions.end()) {
                if (assignments.assigned.size() > incumbent.assigned.size()) {
                    if (proof) {
                        if (params.decide) {
                            vector<pair<NamedVertex, NamedVertex>> solution;
                            for (auto & [l, r] : assignments.assigned)
                                solution.emplace_back(
                                    NamedVertex{l, first.vertex_name(l)},
                                    NamedVertex{r, second.vertex_name(r)});
                            proof->post_solution(solution);
                        }
                        else {
                            vector<tuple<NamedVertex, NamedVertex, bool>> solution;
                            for (int v = 0; v < first.size(); ++v)
                                for (int w = 0; w < second.size(); ++w)
                                    solution.emplace_back(
                                        NamedVertex{v, first.vertex_name(v)},
                                        NamedVertex{w, second.vertex_name(w)},
                                        assignments.assigned.end() != find(assignments.assigned.begin(), assignments.assigned.end(), pair{v, w}));
                            proof->start_level(0);
                            proof->new_incumbent(solution);
                            proof->rewrite_mcs_objective(first.size());
                            proof->start_level(depth);
                        }
                    }

                    if (params.decide) {
                        if (assignments.assigned.size() >= *params.decide) {
                            if (params.count_solutions) {
                                ++solution_count;
                                if (params.enumerate_callback) {
                                    VertexToVertexMapping mapping;
                                    for (auto & [f, s] : assignments.assigned)
                                        mapping.emplace(f, s);
                                    params.enumerate_callback(mapping);
                                }
                                return SearchResult::SatisfiableButKeepGoing;
                            }
                            else {
                                incumbent = assignments;
                                return SearchResult::DecidedTrue;
                            }
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
                    if (proof) {
                        proof->guessing(depth, NamedVertex{left_branch, first.vertex_name(left_branch)},
                            NamedVertex{right_branch, second.vertex_name(right_branch)});
                        proof->start_level(depth + 1);
                    }

                    auto new_domains = branch_assigning(domains, left_branch, right_branch);
                    assignments.assigned.emplace_back(left_branch, right_branch);
                    if (assignments.assigned.size() + bound(new_domains) > incumbent.assigned.size()) {
                        optional<set<int>> new_permitted_branch_variables = nullopt;
                        if (params.connected) {
                            new_permitted_branch_variables = make_optional<set<int>>();
                            if (permitted_branch_variables)
                                new_permitted_branch_variables->insert(permitted_branch_variables->begin(), permitted_branch_variables->end());
                            for (int v = 0; v < first.size(); ++v)
                                if (v != left_branch && first.adjacent(left_branch, v))
                                    new_permitted_branch_variables->emplace(v);
                        }

                        switch (search(depth + 1, assignments, incumbent, new_domains, nodes, solution_count, new_permitted_branch_variables)) {
                        case SearchResult::Aborted: return SearchResult::Aborted;
                        case SearchResult::DecidedTrue: return SearchResult::DecidedTrue;
                        case SearchResult::SatisfiableButKeepGoing: break;
                        case SearchResult::Complete: break;
                        }
                    }
                    else if (proof)
                        proof->mcs_bound(new_domains.partitions);

                    if (proof) {
                        proof->start_level(depth);
                        proof->incorrect_guess(assignments_as_proof_decisions(assignments), true);
                        proof->forget_level(depth + 1);
                    }

                    assignments.assigned.pop_back();
                }

                // now with left_branch assigned to null
                if (proof) {
                    proof->guessing(depth, NamedVertex{left_branch, first.vertex_name(left_branch)},
                        NamedVertex{second.size(), "null"});
                    proof->start_level(depth + 1);
                }

                auto new_domains = branch_rejecting(domains, left_branch);
                assignments.rejected.emplace_back(left_branch);
                if (assignments.assigned.size() + bound(new_domains) > incumbent.assigned.size()) {
                    switch (search(depth + 1, assignments, incumbent, new_domains, nodes, solution_count, permitted_branch_variables)) {
                    case SearchResult::Aborted: return SearchResult::Aborted;
                    case SearchResult::DecidedTrue: return SearchResult::DecidedTrue;
                    case SearchResult::SatisfiableButKeepGoing: break;
                    case SearchResult::Complete: break;
                    }
                }
                else if (proof)
                    proof->mcs_bound(new_domains.partitions);
                if (proof) {
                    proof->start_level(depth);
                    proof->incorrect_guess(assignments_as_proof_decisions(assignments), true);
                    proof->forget_level(depth + 1);
                }
                assignments.rejected.pop_back();
            }

            if (proof)
                proof->out_of_guesses(assignments_as_proof_decisions(assignments));

            return SearchResult::Complete;
        }

        auto run() -> CommonSubgraphResult
        {
            CommonSubgraphResult result;

            map<pair<bool, string_view>, pair<set<int>, set<int>>> initial_partitions;

            for (int v = 0; v < first.size(); ++v)
                initial_partitions[pair{first.adjacent(v, v), first.vertex_label(v)}].first.insert(v);
            for (int v = 0; v < second.size(); ++v)
                initial_partitions[pair{second.adjacent(v, v), second.vertex_label(v)}].second.insert(v);

            SplitDomains domains;
            for (auto & [k, p] : initial_partitions) {
                auto & [l, r] = p;
                if ((! l.empty()) && (! r.empty()))
                    domains.partitions.emplace_back(l, r);
            }

            Assignments assignments, incumbent;

            if (params.decide)
                for (unsigned i = 1; i <= *params.decide - 1; ++i)
                    incumbent.assigned.emplace_back(-int(i), -int(i));

            if (params.decide && (bound(domains) < *params.decide)) {
                result.complete = true;
                if (proof)
                    proof->mcs_bound(domains.partitions);
            }
            else {
                switch (search(0, assignments, incumbent, domains, result.nodes, result.solution_count, nullopt)) {
                case SearchResult::Aborted:
                    break;

                case SearchResult::DecidedTrue:
                    result.complete = true;
                    for (auto & [f, s] : incumbent.assigned)
                        result.mapping.emplace(f, s);
                    break;

                case SearchResult::Complete:
                    result.complete = true;
                    if (! params.decide) {
                        for (auto & [f, s] : incumbent.assigned)
                            result.mapping.emplace(f, s);
                    }
                    break;

                case SearchResult::SatisfiableButKeepGoing:
                    result.complete = true;
                    break;
                }
            }

            if (proof && params.decide && result.complete && result.mapping.empty())
                proof->finish_unsat_proof();
            else if (proof && ! params.decide && result.complete)
                proof->finish_optimisation_proof(first.size() - result.mapping.size());

            return result;
        }
    };
}

auto gss::solve_common_subgraph_problem(const InputGraph & first, const InputGraph & second,
    const CommonSubgraphParams & params) -> CommonSubgraphResult
{
    if (params.count_solutions && ! params.decide)
        throw UnsupportedConfiguration{"Solution counting only makes sense for decision problems"};

    shared_ptr<Proof> proof;
    if (params.proof_options) {
        proof = make_shared<Proof>(*params.proof_options);

        for (int n = 0; n < first.size(); ++n) {
            proof->create_cp_variable(
                n, second.size() + 1,
                [&](int v) { return first.vertex_name(v); },
                [&](int v) { if (v == second.size()) return string("null"); else return second.vertex_name(v); });
        }

        proof->create_null_decision_bound(first.size(), second.size(), params.decide);

        // generate constraints for injectivity
        proof->create_injectivity_constraints(first.size(), second.size());

        // generate edge constraints, and also handle loops here
        for (int p = 0; p < first.size(); ++p) {
            for (int t = 0; t < second.size(); ++t) {
                if (first.adjacent(p, p) != second.adjacent(t, t))
                    proof->create_forbidden_assignment_constraint(p, t);
                else if (first.vertex_label(p) != second.vertex_label(t))
                    proof->create_forbidden_assignment_constraint(p, t);
                else {
                    proof->start_adjacency_constraints_for(p, t);

                    // if p can be mapped to t, then each (non-)neighbour of p...
                    for (int q = 0; q < first.size(); ++q)
                        if (q != p) {
                            // ... must be mapped to a (non-)neighbour of t
                            vector<int> permitted;
                            for (int u = 0; u < second.size(); ++u)
                                if (t != u && first.adjacent(p, q) == second.adjacent(t, u))
                                    permitted.push_back(u);
                            // or null
                            permitted.push_back(second.size());
                            proof->create_adjacency_constraint(p, q, t, permitted, vector<int>{}, false);
                        }
                }
            }
        }

        if (params.connected)
            proof->create_connected_constraints(first.size(), second.size(), [&](int a, int b) { return first.adjacent(a, b); });

        // output the model file
        proof->finalise_model();

        // rewrite the objective to the form we need
        if (! params.decide) {
            vector<tuple<NamedVertex, NamedVertex, bool>> solution;
            for (int v = 0; v < first.size(); ++v)
                for (int w = 0; w < second.size(); ++w)
                    solution.emplace_back(
                        NamedVertex{v, first.vertex_name(v)},
                        NamedVertex{w, second.vertex_name(w)},
                        false);
            proof->new_incumbent(solution);
            proof->rewrite_mcs_objective(first.size());
        }
    }

    if (params.clique) {
        CliqueParams clique_params;
        clique_params.timeout = params.timeout;
        clique_params.start_time = params.start_time;
        clique_params.decide = params.decide;
        clique_params.restarts_schedule = make_unique<NoRestartsSchedule>();
        clique_params.adjust_objective_for_mcs = make_optional(first.size());

        InputGraph assoc{0, false, false};
        vector<pair<int, int>> assoc_encoding, zero_in_proof_objectives;

        for (int v = 0; v < first.size(); ++v)
            for (int w = 0; w < second.size(); ++w)
                if ((first.adjacent(v, v) == second.adjacent(w, w)) && (first.vertex_label(v) == second.vertex_label(w)))
                    assoc_encoding.emplace_back(v, w);
                else if (proof)
                    zero_in_proof_objectives.emplace_back(v, w);

        if (proof)
            proof->create_clique_encoding(assoc_encoding, zero_in_proof_objectives);

        assoc.resize(assoc_encoding.size());
        for (unsigned v = 0; v < assoc_encoding.size(); ++v)
            for (unsigned w = 0; w < assoc_encoding.size(); ++w)
                if (v != w) {
                    auto [vf, vs] = assoc_encoding[v];
                    auto [wf, ws] = assoc_encoding[w];
                    bool edge = false;
                    if (vf != wf && vs != ws && first.adjacent(vf, wf) == second.adjacent(vs, ws)) {
                        if ((! first.adjacent(vf, wf)) || (first.edge_label(vf, wf) == second.edge_label(vs, ws))) {
                            edge = true;
                            assoc.add_edge(v, w);
                        }
                    }

                    if (proof && ! edge)
                        proof->create_clique_nonedge(v, w);
                }

        clique_params.extend_proof = proof;

        if (params.connected) {
            clique_params.connected = [&](int x, const function<auto(int)->int> & invorder) -> SVOBitset {
                auto [f, _] = assoc_encoding[x];
                SVOBitset v(assoc_encoding.size(), 0);
                v.reset();
                for (unsigned y = 0; y < assoc_encoding.size(); ++y)
                    if (first.adjacent(f, assoc_encoding[y].first))
                        v.set(invorder(y));
                return v;
            };
        }

        auto clique_result = solve_clique_problem(assoc, clique_params);

        // now translate the result back into what we expect
        CommonSubgraphResult result;
        for (auto & m : clique_result.clique)
            result.mapping.emplace(assoc_encoding[m]);
        result.nodes = clique_result.nodes;
        result.extra_stats = move(clique_result.extra_stats);
        result.extra_stats.emplace_back("used_clique_solver = true");
        result.complete = clique_result.complete;

        return result;
    }
    else {
        CommonSubgraphRunner runner{first, second, params, proof};
        return runner.run();
    }
}
