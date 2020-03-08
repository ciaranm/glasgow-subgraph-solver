/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "common_subgraph.hh"

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
using std::string_view;
using std::vector;

namespace
{
    enum class SearchResult
    {
        Aborted,
        Complete
    };

    struct SplitDomains
    {
        vector<pair<set<int>, set<int> > > partitions;
    };

    struct Assignments
    {
        vector<pair<int, int> > assigned;
    };

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

        auto bound(const SplitDomains & d) -> int
        {
            int result = 0;
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
                Assignments & assignments,
                Assignments & incumbent,
                const SplitDomains & domains,
                unsigned long long & nodes) -> SearchResult
        {
            ++nodes;

            auto branch_domains = find_branch_partition(domains);
            if (branch_domains == domains.partitions.end()) {
                if (assignments.assigned.size() > incumbent.assigned.size())
                    incumbent = assignments;
            }
            else {
                int left_branch = *branch_domains->first.begin();
                for (auto & right_branch : branch_domains->second) {
                    // branch with left_branch assigned to right_branch
                    auto new_domains = branch_assigning(domains, left_branch, right_branch);
                    if (assignments.assigned.size() + bound(new_domains) + 1 > incumbent.assigned.size()) {
                        assignments.assigned.emplace_back(left_branch, right_branch);
                        if (SearchResult::Aborted == search(assignments, incumbent, new_domains, nodes))
                            return SearchResult::Aborted;
                        assignments.assigned.pop_back();
                    }
                }

                // now with left_branch assigned to null
                auto new_domains = branch_rejecting(domains, left_branch);
                if (assignments.assigned.size() + bound(new_domains) > incumbent.assigned.size()) {
                    if (SearchResult::Aborted == search(assignments, incumbent, new_domains, nodes))
                        return SearchResult::Aborted;
                }
            }

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
            switch (search(assignments, incumbent, domains, result.nodes)) {
                case SearchResult::Aborted:
                    break;

                case SearchResult::Complete:
                    result.complete = true;
                    for (auto & [ f, s ] : incumbent.assigned)
                        result.mapping.emplace(f, s);
                    break;
            }

            return result;
        }
    };
}

auto solve_common_subgraph_problem(const InputGraph & first, const InputGraph & second, const CommonSubgraphParams & params) -> CommonSubgraphResult
{
    CommonSubgraphRunner runner{ first, second, params };
    return runner.run();
}

