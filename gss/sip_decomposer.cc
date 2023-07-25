#include <gss/homomorphism.hh>
#include <gss/innards/proof.hh>
#include <gss/loooong.hh>
#include <gss/sip_decomposer.hh>

#include <numeric>
#include <set>
#include <string>
#include <vector>

using namespace gss;
using namespace gss::innards;

using std::gcd;
using std::pair;
using std::set;
using std::to_string;
using std::vector;

namespace
{
    auto find_removable_isolated_pattern_vertices(
        const InputGraph & pattern,
        const HomomorphismParams & params,
        set<int> & isolated_pattern_vertices) -> void
    {
        if (params.induced || (! params.pattern_less_constraints.empty()) || (! params.target_occur_less_constraints.empty()) || params.lackey || pattern.has_vertex_labels() || (params.count_solutions && params.enumerate_callback) || pattern.directed())
            return;

        for (int i = 0; i < pattern.size(); ++i)
            if (0 == pattern.degree(i))
                isolated_pattern_vertices.emplace(i);
    }

    template <typename I_>
    auto n_choose_k(I_ n, I_ k) -> I_
    {
        I_ r = 1;
        for (I_ d = 1; d <= k; ++d) {
            I_ numerator = n - k + d;
            I_ denominator = d;
            I_ cf = gcd(r, denominator);
            r /= cf;
            denominator /= cf;
            numerator /= denominator;
            r *= numerator;
        }
        return r;
    }

    template <typename I_>
    auto factorial(I_ n) -> I_
    {
        I_ r = 1;
        for (I_ d = 2; d <= n; ++d)
            r *= d;
        return r;
    }
}

auto gss::solve_sip_by_decomposition(const InputGraph & pattern, const InputGraph & target,
    const HomomorphismParams & params) -> HomomorphismResult
{
    set<int> isolated_pattern_vertices;
    find_removable_isolated_pattern_vertices(pattern, params, isolated_pattern_vertices);

    if (! isolated_pattern_vertices.empty()) {
        InputGraph reduced_pattern(pattern.size() - isolated_pattern_vertices.size(), pattern.has_vertex_labels(),
            pattern.has_edge_labels());

        vector<int> original_to_reduced, reduced_to_original;

        for (int i = 0; i < pattern.size(); ++i) {
            if (isolated_pattern_vertices.count(i))
                original_to_reduced.push_back(-1);
            else {
                original_to_reduced.push_back(reduced_to_original.size());
                reduced_to_original.push_back(i);
            }
        }

        for (int i = 0; i < pattern.size(); ++i) {
            auto r_i = original_to_reduced.at(i);
            if (r_i == -1)
                continue;

            for (int j = 0; j < pattern.size(); ++j) {
                auto r_j = original_to_reduced.at(j);
                if (r_j == -1)
                    continue;

                if (pattern.adjacent(i, j))
                    reduced_pattern.add_directed_edge(r_i, r_j, "");
            }
        }
        auto result = solve_homomorphism_problem(reduced_pattern, target, params);

        result.extra_stats.emplace_back("isolated_pattern_vertices = " + to_string(isolated_pattern_vertices.size()));

        // fix up the result to be in the non-decomposed form
        if (! result.mapping.empty()) {
            // need to re-add isolated vertices, arbitrarily. what's available?
            set<int> t_avail;
            for (int i = 0; i < target.size(); ++i)
                t_avail.emplace(i);

            // convert back to original indices, and remove used vertices
            auto sub_mapping = result.mapping;
            result.mapping.clear();
            for (auto & [p, t] : sub_mapping) {
                t_avail.erase(t);
                result.mapping.emplace(reduced_to_original.at(p), t);
            }

            // add in isolated vertices
            for (auto & p : isolated_pattern_vertices) {
                result.mapping.emplace(p, *t_avail.begin());
                t_avail.erase(t_avail.begin());
            }
        }

        // fix up the solution count
        if (params.count_solutions && result.solution_count > 0) {
            loooong unmapped_target_vertices = target.size() - reduced_pattern.size();
            loooong solution_multiplier = n_choose_k<loooong>(unmapped_target_vertices, isolated_pattern_vertices.size());
            loooong isolated_symmetry_multiplier = factorial<loooong>(isolated_pattern_vertices.size());
            result.solution_count *= solution_multiplier * isolated_symmetry_multiplier;
        }

        return result;
    }
    else
        return solve_homomorphism_problem(pattern, target, params);
}
