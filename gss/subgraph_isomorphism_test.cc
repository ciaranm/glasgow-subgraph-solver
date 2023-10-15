#include <gss/formats/csv.hh>
#include <gss/homomorphism.hh>

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
