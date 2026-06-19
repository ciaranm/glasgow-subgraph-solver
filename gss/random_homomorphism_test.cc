#include <gss/formats/input_graph.hh>
#include <gss/homomorphism.hh>
#include <gss/innards/verify.hh>

#include <catch2/catch_test_macros.hpp>

#include <chrono>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <vector>

using namespace gss;
using namespace gss::innards;

using std::make_shared;
using std::make_unique;
using std::map;
using std::mt19937;
using std::set;
using std::uniform_real_distribution;
using std::vector;

using std::chrono::operator""s;

namespace
{
    // A random undirected graph with optional self-loops. (add_edge is symmetric, so
    // there is no directed-graph API to drive here; directed instances are left for a
    // future extension.)
    auto random_graph(int n, double edge_probability, double loop_probability, mt19937 & rng) -> InputGraph
    {
        InputGraph g{n, false, false};
        uniform_real_distribution<double> dist{0.0, 1.0};
        for (int v = 0; v < n; ++v) {
            if (dist(rng) < loop_probability)
                g.add_edge(v, v);
            for (int w = v + 1; w < n; ++w)
                if (dist(rng) < edge_probability)
                    g.add_edge(v, w);
        }
        return g;
    }

    // The brute-force oracle: every complete mapping pattern -> target that the
    // solver's own verifier accepts. verify_homomorphism does no search, so this is
    // an independent specification of "what counts as a solution" for these options.
    auto brute_force_solutions(const InputGraph & pattern, const InputGraph & target,
        bool injective, bool locally_injective, bool induced) -> set<map<int, int>>
    {
        set<map<int, int>> solutions;
        int np = pattern.size(), nt = target.size();
        if (np == 0 || nt == 0)
            return solutions;

        // iterate over all nt^np functions as a mixed-radix counter
        vector<int> assignment(np, 0);
        for (;;) {
            map<int, int> mapping;
            for (int i = 0; i < np; ++i)
                mapping.emplace(i, assignment[i]);

            try {
                verify_homomorphism(pattern, target, injective, locally_injective, induced, mapping);
                solutions.insert(move(mapping));
            }
            catch (const BuggySolution &) {
                // not a valid mapping for these options
            }

            int i = 0;
            for (; i < np; ++i) {
                if (++assignment[i] < nt)
                    break;
                assignment[i] = 0;
            }
            if (i == np)
                break;
        }

        return solutions;
    }
}

TEST_CASE("random instances: solver enumeration matches the brute-force oracle")
{
    // Fixed seed so failures reproduce; the iteration index identifies the instance.
    mt19937 rng{0x5eed};

    for (int iter = 0; iter < 1500; ++iter) {
        int np = 1 + int(rng() % 4); // 1..4
        int nt = 1 + int(rng() % 5); // 1..5
        double pattern_density = (rng() % 100) / 100.0;
        double target_density = (rng() % 100) / 100.0;
        double loop_probability = (rng() % 3 == 0) ? 0.25 : 0.0;

        auto pattern = random_graph(np, pattern_density, loop_probability, rng);
        auto target = random_graph(nt, target_density, loop_probability, rng);

        bool induced = (rng() % 2);
        int injectivity_choice = rng() % 3;
        auto injectivity = injectivity_choice == 0 ? Injectivity::Injective
            : injectivity_choice == 1 ? Injectivity::NonInjective
                                      : Injectivity::LocallyInjective;
        bool injective = (injectivity == Injectivity::Injective);
        bool locally_injective = (injectivity == Injectivity::LocallyInjective);

        auto expected = brute_force_solutions(pattern, target, injective, locally_injective, induced);

        HomomorphismParams params;
        params.timeout = make_shared<Timeout>(0s);
        params.restarts_schedule = make_unique<NoRestartsSchedule>();
        params.injectivity = injectivity;
        params.induced = induced;
        params.count_solutions = true;

        set<map<int, int>> got;
        params.enumerate_callback = [&](const VertexToVertexMapping & mapping) {
            got.insert(mapping);
            return true;
        };

        // TEMPORARY DIAGNOSTIC (debug/li-sigill-24.04): flush the instance descriptor
        // before solving, so a fatal signal in solve identifies the exact instance.
        {
            auto dbg_edges = [](const InputGraph & g) {
                std::string s;
                for (int a = 0; a < g.size(); ++a)
                    for (int b = 0; b < g.size(); ++b)
                        if (g.adjacent(a, b))
                            s += " " + std::to_string(a) + "-" + std::to_string(b);
                return s;
            };
            std::cerr << "ITER " << iter << " np=" << np << " nt=" << nt
                      << " induced=" << induced << " inj=" << injectivity_choice
                      << " P[" << dbg_edges(pattern) << " ] T[" << dbg_edges(target) << " ]"
                      << std::endl;
        }

        auto result = solve_homomorphism_problem(pattern, target, params);

        auto edges = [](const InputGraph & g) {
            std::string s;
            for (int a = 0; a < g.size(); ++a)
                for (int b = a; b < g.size(); ++b)
                    if (g.adjacent(a, b))
                        s += " " + std::to_string(a) + "-" + std::to_string(b);
            return s;
        };
        INFO("iteration " << iter << ", |pattern|=" << np << ", |target|=" << nt
                          << ", induced=" << induced << ", injectivity=" << injectivity_choice
                          << "\n  pattern edges:" << edges(pattern)
                          << "\n  target edges:" << edges(target)
                          << "\n  solver=" << got.size() << " oracle=" << expected.size());
        CHECK(got == expected);
        CHECK(result.solution_count == int(expected.size()));
        // a non-empty enumeration must have run to completion; the only incomplete
        // case is the trivial injective pattern-bigger-than-target shortcut (0 solutions)
        CHECK((result.complete || expected.empty()));
    }
}
