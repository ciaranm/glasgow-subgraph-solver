#include <gss/clique.hh>
#include <gss/formats/input_graph.hh>
#include <gss/restarts.hh>
#include <gss/timeout.hh>

#include <catch2/catch_test_macros.hpp>

#include <chrono>
#include <memory>
#include <set>

using namespace gss;

using std::make_shared;
using std::make_unique;
using std::set;

using std::chrono::operator""s;
using std::chrono::steady_clock;

namespace
{
    auto make_params() -> CliqueParams
    {
        CliqueParams params;
        params.timeout = make_shared<Timeout>(0s);
        params.start_time = steady_clock::now();
        params.restarts_schedule = make_unique<NoRestartsSchedule>();
        return params;
    }

    auto is_clique(const InputGraph & g, const set<int> & vs) -> bool
    {
        for (auto a : vs)
            for (auto b : vs)
                if (a != b && ! g.adjacent(a, b))
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
}

TEST_CASE("maximum clique of a complete graph is the whole graph")
{
    auto g = complete_graph(4);
    auto result = solve_clique_problem(g, make_params());
    CHECK(result.clique == set<int>{0, 1, 2, 3});
    CHECK(is_clique(g, result.clique));
}

TEST_CASE("maximum clique of a triangle")
{
    auto g = complete_graph(3);
    auto result = solve_clique_problem(g, make_params());
    CHECK(result.clique.size() == 3);
    CHECK(is_clique(g, result.clique));
}

TEST_CASE("maximum clique of a path is a single edge")
{
    auto g = path_graph(4);
    auto result = solve_clique_problem(g, make_params());
    CHECK(result.clique.size() == 2);
    CHECK(is_clique(g, result.clique));
}

TEST_CASE("maximum clique of an edgeless graph is a single vertex")
{
    InputGraph g{3, false, false};
    auto result = solve_clique_problem(g, make_params());
    CHECK(result.clique.size() == 1);
}

TEST_CASE("maximum clique of two disjoint triangles is one triangle")
{
    InputGraph g{6, false, false};
    for (int base : {0, 3}) {
        g.add_edge(base, base + 1);
        g.add_edge(base + 1, base + 2);
        g.add_edge(base, base + 2);
    }
    auto result = solve_clique_problem(g, make_params());
    CHECK(result.clique.size() == 3);
    CHECK(is_clique(g, result.clique));
}

TEST_CASE("a unique maximum clique is found exactly")
{
    // K4 on {0,1,2,3}, plus vertex 4 attached only to vertex 0. The only
    // 4-clique is {0,1,2,3}; vertex 4 can only ever be in a 2-clique.
    auto g = complete_graph(4);
    g.resize(5);
    g.add_edge(0, 4);
    auto result = solve_clique_problem(g, make_params());
    CHECK(result.clique == set<int>{0, 1, 2, 3});
}

TEST_CASE("decision: a clique of the requested size exists")
{
    auto g = complete_graph(3);
    auto params = make_params();
    params.decide = 3;
    auto result = solve_clique_problem(g, params);
    CHECK(result.clique.size() == 3);
    CHECK(is_clique(g, result.clique));
}

TEST_CASE("decision: a clique of the requested size does not exist")
{
    auto g = complete_graph(3); // no 4-clique
    auto params = make_params();
    params.decide = 4;
    auto result = solve_clique_problem(g, params);
    CHECK(result.clique.empty());
}

TEST_CASE("stop_after_finding returns a clique at least that big")
{
    auto g = complete_graph(5);
    auto params = make_params();
    params.stop_after_finding = 3;
    auto result = solve_clique_problem(g, params);
    CHECK(result.clique.size() >= 3);
    CHECK(is_clique(g, result.clique));
}

TEST_CASE("the colour orderings and vertex order all find the same maximum")
{
    // K4 on {0,1,2,3} plus a pendant: unique maximum clique size is 4.
    auto g = complete_graph(4);
    g.resize(5);
    g.add_edge(0, 4);

    for (auto order : {ColourClassOrder::ColourOrder, ColourClassOrder::SingletonsFirst, ColourClassOrder::Sorted}) {
        for (bool input_order : {false, true}) {
            auto params = make_params();
            params.colour_class_order = order;
            params.input_order = input_order;
            auto result = solve_clique_problem(g, params);
            CHECK(result.clique.size() == 4);
            CHECK(is_clique(g, result.clique));
        }
    }
}
