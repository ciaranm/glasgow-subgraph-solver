#include <gss/common_subgraph.hh>
#include <gss/configuration.hh>
#include <gss/formats/input_graph.hh>
#include <gss/timeout.hh>

#include <catch2/catch_test_macros.hpp>

#include <chrono>
#include <memory>

using namespace gss;

using std::make_shared;

using std::chrono::operator""s;
using std::chrono::steady_clock;

namespace
{
    auto make_params() -> CommonSubgraphParams
    {
        CommonSubgraphParams params;
        params.timeout = make_shared<Timeout>(0s);
        params.start_time = steady_clock::now();
        return params;
    }

    // A mapping is a common *induced* subgraph iff it preserves both edges and
    // non-edges between every pair of mapped vertices.
    auto is_common_induced(const InputGraph & f, const InputGraph & s, const VertexToVertexMapping & m) -> bool
    {
        for (auto & [a1, b1] : m)
            for (auto & [a2, b2] : m)
                if (a1 != a2 && f.adjacent(a1, a2) != s.adjacent(b1, b2))
                    return false;
        return true;
    }

    auto complete_graph(int n) -> InputGraph
    {
        InputGraph g{n, false, false};
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                g.add_edge(i, j);
        return g;
    }

    auto path_graph(int n) -> InputGraph
    {
        InputGraph g{n, false, false};
        for (int i = 0; i + 1 < n; ++i)
            g.add_edge(i, i + 1);
        return g;
    }

    auto two_disjoint_edges() -> InputGraph
    {
        InputGraph g{4, false, false};
        g.add_edge(0, 1);
        g.add_edge(2, 3);
        return g;
    }
}

TEST_CASE("maximum common subgraph of two identical triangles is the whole triangle")
{
    auto f = complete_graph(3), s = complete_graph(3);
    auto result = solve_common_subgraph_problem(f, s, make_params());
    CHECK(result.mapping.size() == 3);
    CHECK(is_common_induced(f, s, result.mapping));
}

TEST_CASE("a triangle is a common subgraph of a larger complete graph")
{
    auto f = complete_graph(3), s = complete_graph(4);
    auto result = solve_common_subgraph_problem(f, s, make_params());
    CHECK(result.mapping.size() == 3);
    CHECK(is_common_induced(f, s, result.mapping));
}

TEST_CASE("induced: a triangle and a path share only an edge")
{
    // K3 has all three edges; P3 has a non-edge, so the largest common *induced*
    // subgraph is a single edge (two vertices).
    auto f = complete_graph(3), s = path_graph(3);
    auto result = solve_common_subgraph_problem(f, s, make_params());
    CHECK(result.mapping.size() == 2);
    CHECK(is_common_induced(f, s, result.mapping));
}

TEST_CASE("two identical paths share the whole path")
{
    auto f = path_graph(4), s = path_graph(4);
    auto result = solve_common_subgraph_problem(f, s, make_params());
    CHECK(result.mapping.size() == 4);
    CHECK(is_common_induced(f, s, result.mapping));
}

TEST_CASE("edgeless graphs share min(|V|) vertices")
{
    InputGraph f{3, false, false}, s{4, false, false};
    auto result = solve_common_subgraph_problem(f, s, make_params());
    CHECK(result.mapping.size() == 3);
    CHECK(is_common_induced(f, s, result.mapping));
}

TEST_CASE("decision: a common subgraph of the requested size exists")
{
    auto f = complete_graph(3), s = complete_graph(3);
    auto params = make_params();
    params.decide = 3;
    auto result = solve_common_subgraph_problem(f, s, params);
    CHECK(result.mapping.size() == 3);
    CHECK(is_common_induced(f, s, result.mapping));
}

TEST_CASE("decision: a common subgraph of the requested size does not exist")
{
    auto f = complete_graph(3), s = complete_graph(3); // no common subgraph of size 4
    auto params = make_params();
    params.decide = 4;
    auto result = solve_common_subgraph_problem(f, s, params);
    CHECK(result.mapping.empty());
}

TEST_CASE("the clique reduction finds the same maximum")
{
    auto f = complete_graph(3), s = complete_graph(4);
    auto params = make_params();
    params.clique = true;
    auto result = solve_common_subgraph_problem(f, s, params);
    CHECK(result.mapping.size() == 3);
    CHECK(is_common_induced(f, s, result.mapping));
}

TEST_CASE("connected vs unconnected common subgraph")
{
    auto f = two_disjoint_edges(), s = two_disjoint_edges();

    SECTION("unconnected: both edges can be mapped")
    {
        auto result = solve_common_subgraph_problem(f, s, make_params());
        CHECK(result.mapping.size() == 4);
        CHECK(is_common_induced(f, s, result.mapping));
    }

    SECTION("connected: only a single edge")
    {
        auto params = make_params();
        params.connected = true;
        auto result = solve_common_subgraph_problem(f, s, params);
        CHECK(result.mapping.size() == 2);
        CHECK(is_common_induced(f, s, result.mapping));
    }
}

TEST_CASE("counting solutions only makes sense as a decision problem")
{
    auto f = complete_graph(3), s = complete_graph(3);
    auto params = make_params();
    params.count_solutions = true; // but no params.decide
    CHECK_THROWS_AS(solve_common_subgraph_problem(f, s, params), UnsupportedConfiguration);
}

TEST_CASE("counting common subgraphs of a given size")
{
    // K2 -> K2: the two size-2 mappings are the two bijections of an edge.
    auto f = complete_graph(2), s = complete_graph(2);
    auto params = make_params();
    params.decide = 2;
    params.count_solutions = true;
    auto result = solve_common_subgraph_problem(f, s, params);
    CHECK(result.solution_count == 2);
}
