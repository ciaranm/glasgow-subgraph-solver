#include <gss/formats/csv.hh>
#include <gss/homomorphism.hh>
#include <gss/loooong.hh>
#include <gss/sip_decomposer.hh>

#include <catch2/catch_test_macros.hpp>

#include <sstream>
#include <utility>

using namespace gss;

using std::chrono::operator""s;
using std::make_shared;
using std::make_unique;
using std::stringstream;

TEST_CASE("subgraph isomorphism no edges")
{
    auto pattern = read_csv(stringstream{// clang-format off
R"(a,
b,
c,
)"}, "pattern"); // clang-format on

    auto target = read_csv(stringstream{// clang-format off
R"(1,
2,
3,
4,
)"}, "target"); // clang-format on

    HomomorphismParams params;
    params.timeout = make_shared<Timeout>(0s);
    params.restarts_schedule = make_unique<NoRestartsSchedule>();

    SECTION("decide")
    {
        auto result = solve_homomorphism_problem(pattern, target, params);
        CHECK(result.mapping.size() == 3);
        CHECK(result.complete);
    }

    SECTION("count")
    {
        params.count_solutions = true;
        auto result = solve_homomorphism_problem(pattern, target, params);
        CHECK(result.solution_count == 4 * 3 * 2);
        CHECK(result.complete);
    }
}

TEST_CASE("subgraph isomorphism line")
{
    auto pattern = read_csv(stringstream{// clang-format off
R"(a,b
b,c
)"}, "pattern"); // clang-format on

    auto target = read_csv(stringstream{// clang-format off
R"(1,2
2,3
3,4
)"}, "target"); // clang-format on

    HomomorphismParams params;
    params.timeout = make_shared<Timeout>(0s);
    params.restarts_schedule = make_unique<NoRestartsSchedule>();

    SECTION("count")
    {
        params.count_solutions = true;
        auto result = solve_homomorphism_problem(pattern, target, params);
        CHECK(result.solution_count == 4);
        CHECK(result.complete);
    }
}

TEST_CASE("subgraph isomorphism loop")
{
    auto pattern = read_csv(stringstream{// clang-format off
R"(a,a
a,b
)"}, "pattern"); // clang-format on

    auto target = read_csv(stringstream{// clang-format off
R"(1,2
2,2
2,3
3,4
)"}, "target"); // clang-format on

    HomomorphismParams params;
    params.timeout = make_shared<Timeout>(0s);
    params.restarts_schedule = make_unique<NoRestartsSchedule>();

    SECTION("count")
    {
        params.count_solutions = true;
        auto result = solve_homomorphism_problem(pattern, target, params);
        CHECK(result.solution_count == 2);
        CHECK(result.complete);
    }
}

TEST_CASE("subgraph isomorphism with undirected edge labels")
{
    // A two-edge path, one red edge and one blue.
    auto pattern = read_csv(stringstream{// clang-format off
R"(a,b,red
b,c,blue
)"}, "pattern"); // clang-format on

    // A four-vertex path: red, blue, red.
    auto target = read_csv(stringstream{// clang-format off
R"(1,2,red
2,3,blue
3,4,red
)"}, "target"); // clang-format on

    CHECK(pattern.has_edge_labels());
    CHECK_FALSE(pattern.directed()); // labelling an undirected edge must not make it directed

    HomomorphismParams params;
    params.timeout = make_shared<Timeout>(0s);
    params.restarts_schedule = make_unique<NoRestartsSchedule>();
    params.count_solutions = true;

    // b is the only pattern vertex on both a red and a blue edge, so it maps to 2 or
    // to 3, and either choice fixes the rest: (a, b, c) -> (1, 2, 3) or (4, 3, 2).
    SECTION("count")
    {
        auto result = solve_homomorphism_problem(pattern, target, params);
        CHECK(result.solution_count == 2);
        CHECK(result.complete);
    }

    // a and c are non-adjacent, and so are their images in both mappings above.
    SECTION("induced")
    {
        params.induced = true;
        auto result = solve_homomorphism_problem(pattern, target, params);
        CHECK(result.solution_count == 2);
        CHECK(result.complete);
    }
}

TEST_CASE("subgraph isomorphism with a uniform edge label matches the unlabelled count")
{
    // The graphs from "subgraph isomorphism line" with every edge labelled the same.
    // The labels rule nothing out, so the count must not move.
    auto pattern = read_csv(stringstream{// clang-format off
R"(a,b,red
b,c,red
)"}, "pattern"); // clang-format on

    auto target = read_csv(stringstream{// clang-format off
R"(1,2,red
2,3,red
3,4,red
)"}, "target"); // clang-format on

    HomomorphismParams params;
    params.timeout = make_shared<Timeout>(0s);
    params.restarts_schedule = make_unique<NoRestartsSchedule>();
    params.count_solutions = true;

    auto result = solve_homomorphism_problem(pattern, target, params);
    CHECK(result.solution_count == 4);
    CHECK(result.complete);
}

TEST_CASE("decomposition keeps the edge labels of the reduced pattern")
{
    // One labelled edge plus an isolated vertex, so the decomposition actually
    // triggers. Rebuilding the reduced pattern must not lose the label, or the
    // reduced problem asks for empty-labelled target edges and finds nothing.
    auto pattern = read_csv(stringstream{// clang-format off
R"(a,b,red
z,
)"}, "pattern"); // clang-format on

    auto target = read_csv(stringstream{// clang-format off
R"(1,2,red
2,3,red
9,
)"}, "target"); // clang-format on

    HomomorphismParams params;
    params.timeout = make_shared<Timeout>(0s);
    params.restarts_schedule = make_unique<NoRestartsSchedule>();
    params.count_solutions = true;

    auto undecomposed = solve_homomorphism_problem(pattern, target, params);
    REQUIRE(undecomposed.complete);
    CHECK(undecomposed.solution_count == 8);

    auto decomposed = solve_sip_by_decomposition(pattern, target, params);
    CHECK(decomposed.complete);
    CHECK(decomposed.solution_count == undecomposed.solution_count);
}
