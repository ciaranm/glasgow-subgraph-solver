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

    // A mapping is a common *induced* subgraph iff, for every ordered pair of mapped
    // vertices, adjacency agrees on both sides. Taking the pair (a, a) into account
    // (rather than skipping it) also enforces the loop rule: a vertex may be mapped
    // to another only if either both have a loop or neither does.
    auto is_common_induced(const InputGraph & f, const InputGraph & s, const VertexToVertexMapping & m) -> bool
    {
        for (auto & [a1, b1] : m)
            for (auto & [a2, b2] : m)
                if (f.adjacent(a1, a2) != s.adjacent(b1, b2))
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

// ---------------------------------------------------------------------------
// Loops
// ---------------------------------------------------------------------------

TEST_CASE("a looped vertex may only be associated with a looped vertex")
{
    InputGraph looped{1, false, false};
    looped.add_edge(0, 0);
    InputGraph plain{1, false, false};

    auto result = solve_common_subgraph_problem(looped, plain, make_params());
    CHECK(result.mapping.empty()); // loop vs non-loop -> nothing in common
}

TEST_CASE("the maximum common subgraph respects loops")
{
    // first: edge 0-1 with a loop on 0; second: edge 0-1 with a loop on 1.
    InputGraph f{2, false, false};
    f.add_edge(0, 1);
    f.add_edge(0, 0);
    InputGraph s{2, false, false};
    s.add_edge(0, 1);
    s.add_edge(1, 1);

    auto result = solve_common_subgraph_problem(f, s, make_params());
    CHECK(result.mapping.size() == 2);
    CHECK(is_common_induced(f, s, result.mapping));
    // the only loop-respecting mapping pairs the two loopy vertices.
    CHECK(result.mapping.at(0) == 1);
    CHECK(result.mapping.at(1) == 0);
}

TEST_CASE("is_common_induced rejects a loop-rule violation")
{
    // Guards the checker itself: mapping a looped vertex onto a non-looped one
    // is not a common induced subgraph.
    InputGraph looped{1, false, false};
    looped.add_edge(0, 0);
    InputGraph plain{1, false, false};
    CHECK_FALSE(is_common_induced(looped, plain, VertexToVertexMapping{{0, 0}}));
}

// ---------------------------------------------------------------------------
// Labels
// ---------------------------------------------------------------------------

TEST_CASE("vertex labels constrain the association")
{
    SECTION("matching labels allow a full mapping")
    {
        InputGraph f{2, true, false};
        f.set_vertex_label(0, "x");
        f.set_vertex_label(1, "y");
        InputGraph s{2, true, false};
        s.set_vertex_label(0, "x");
        s.set_vertex_label(1, "y");

        auto result = solve_common_subgraph_problem(f, s, make_params());
        CHECK(result.mapping.size() == 2);
        CHECK(f.vertex_label(0) == s.vertex_label(result.mapping.at(0)));
    }

    SECTION("a missing label reduces the mapping")
    {
        InputGraph f{2, true, false};
        f.set_vertex_label(0, "x");
        f.set_vertex_label(1, "y");
        InputGraph s{2, true, false};
        s.set_vertex_label(0, "x");
        s.set_vertex_label(1, "x"); // no "y" available

        auto result = solve_common_subgraph_problem(f, s, make_params());
        CHECK(result.mapping.size() == 1); // only the "x" vertex can be placed
    }
}

TEST_CASE("edge labels constrain the association")
{
    SECTION("matching edge labels allow the edge to be mapped")
    {
        InputGraph f{2, false, true};
        f.add_directed_edge(0, 1, "a");
        f.add_directed_edge(1, 0, "a");
        InputGraph s{2, false, true};
        s.add_directed_edge(0, 1, "a");
        s.add_directed_edge(1, 0, "a");

        auto result = solve_common_subgraph_problem(f, s, make_params());
        CHECK(result.mapping.size() == 2);
    }

    SECTION("differing edge labels prevent the edge being mapped")
    {
        InputGraph f{2, false, true};
        f.add_directed_edge(0, 1, "a");
        f.add_directed_edge(1, 0, "a");
        InputGraph s{2, false, true};
        s.add_directed_edge(0, 1, "b");
        s.add_directed_edge(1, 0, "b");

        auto result = solve_common_subgraph_problem(f, s, make_params());
        CHECK(result.mapping.size() == 1); // the "a" edge has no "a" edge to map to
    }
}

TEST_CASE("the clique reduction also respects vertex labels")
{
    InputGraph f{2, true, false};
    f.set_vertex_label(0, "x");
    f.set_vertex_label(1, "y");
    InputGraph s{2, true, false};
    s.set_vertex_label(0, "x");
    s.set_vertex_label(1, "x");

    auto params = make_params();
    params.clique = true;
    auto result = solve_common_subgraph_problem(f, s, params);
    CHECK(result.mapping.size() == 1); // same as the non-clique path
}
