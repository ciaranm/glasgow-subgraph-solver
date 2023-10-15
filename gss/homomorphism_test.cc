#include <gss/formats/csv.hh>
#include <gss/homomorphism.hh>

#include <catch2/catch_test_macros.hpp>

#include <set>
#include <sstream>
#include <utility>

using namespace gss;

using std::chrono::operator""s;
using std::make_shared;
using std::make_unique;
using std::pair;
using std::set;
using std::string;
using std::stringstream;

TEST_CASE("homomorphism no edges")
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
    params.injectivity = Injectivity::NonInjective;

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
        CHECK(result.solution_count == 4 * 4 * 4);
        CHECK(result.complete);
    }
}

TEST_CASE("homomorphism line")
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
    params.injectivity = Injectivity::NonInjective;

    auto p = [&](const string & s) { return *pattern.vertex_from_name(s); };
    auto t = [&](const string & s) { return *target.vertex_from_name(s); };

    SECTION("count")
    {
        params.count_solutions = true;
        auto result = solve_homomorphism_problem(pattern, target, params);
        CHECK(result.solution_count == 10);
        CHECK(result.complete);
    }

    SECTION("enumerate")
    {
        set<VertexToVertexMapping> got;
        params.count_solutions = true;
        params.enumerate_callback = [&](const VertexToVertexMapping & m) -> bool {
            got.insert(m);
            return true;
        };
        auto result = solve_homomorphism_problem(pattern, target, params);
        CHECK(got == set<VertexToVertexMapping>{
                         {{p("a"), t("1")}, {p("b"), t("2")}, {p("c"), t("1")}}, //
                         {{p("a"), t("1")}, {p("b"), t("2")}, {p("c"), t("3")}}, //
                         {{p("a"), t("2")}, {p("b"), t("1")}, {p("c"), t("2")}}, //
                         {{p("a"), t("2")}, {p("b"), t("3")}, {p("c"), t("2")}}, //
                         {{p("a"), t("2")}, {p("b"), t("3")}, {p("c"), t("4")}}, //
                         {{p("a"), t("3")}, {p("b"), t("2")}, {p("c"), t("1")}}, //
                         {{p("a"), t("3")}, {p("b"), t("2")}, {p("c"), t("3")}}, //
                         {{p("a"), t("3")}, {p("b"), t("4")}, {p("c"), t("3")}}, //
                         {{p("a"), t("4")}, {p("b"), t("3")}, {p("c"), t("2")}}, //
                         {{p("a"), t("4")}, {p("b"), t("3")}, {p("c"), t("4")}}  //
                     });
    }
}

TEST_CASE("homomorphism loop")
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
    params.injectivity = Injectivity::NonInjective;

    auto p = [&](const string & s) { return *pattern.vertex_from_name(s); };
    auto t = [&](const string & s) { return *target.vertex_from_name(s); };

    SECTION("decide")
    {
        auto result = solve_homomorphism_problem(pattern, target, params);
        CHECK(result.solution_count == 1);
        CHECK(result.nodes == 0);
        CHECK(result.complete);
    }

    SECTION("count")
    {
        params.count_solutions = true;
        auto result = solve_homomorphism_problem(pattern, target, params);
        CHECK(result.solution_count == 3);
        CHECK(result.complete);
    }

    SECTION("enumerate")
    {
        set<VertexToVertexMapping> got;
        params.count_solutions = true;
        params.enumerate_callback = [&](const VertexToVertexMapping & m) -> bool {
            got.insert(m);
            return true;
        };
        auto result = solve_homomorphism_problem(pattern, target, params);
        CHECK(got == set<VertexToVertexMapping>{
                         {{p("a"), t("2")}, {p("b"), t("1")}}, //
                         {{p("a"), t("2")}, {p("b"), t("2")}}, //
                         {{p("a"), t("2")}, {p("b"), t("3")}}  //
                     });
    }
}

TEST_CASE("homomorphism loop shrinking")
{
    auto pattern = read_csv(stringstream{// clang-format off
R"(a,a
a,b
b,c
)"}, "pattern"); // clang-format on

    auto target = read_csv(stringstream{// clang-format off
R"(1,1
1,2
)"}, "target"); // clang-format on

    HomomorphismParams params;
    params.timeout = make_shared<Timeout>(0s);
    params.restarts_schedule = make_unique<NoRestartsSchedule>();
    params.injectivity = Injectivity::NonInjective;

    auto p = [&](const string & s) { return *pattern.vertex_from_name(s); };
    auto t = [&](const string & s) { return *target.vertex_from_name(s); };

    SECTION("decide")
    {
        auto result = solve_homomorphism_problem(pattern, target, params);
        CHECK(result.solution_count == 1);
        CHECK(result.nodes == 0);
        CHECK(result.complete);
    }

    SECTION("count")
    {
        params.count_solutions = true;
        auto result = solve_homomorphism_problem(pattern, target, params);
        CHECK(result.solution_count == 3);
        CHECK(result.complete);
    }

    SECTION("enumerate")
    {
        set<VertexToVertexMapping> got;
        params.count_solutions = true;
        params.enumerate_callback = [&](const VertexToVertexMapping & m) -> bool {
            got.insert(m);
            return true;
        };
        auto result = solve_homomorphism_problem(pattern, target, params);
        CHECK(got == set<VertexToVertexMapping>{
                         {{p("a"), t("1")}, {p("b"), t("1")}, {p("c"), t("1")}}, //
                         {{p("a"), t("1")}, {p("b"), t("1")}, {p("c"), t("2")}}, //
                         {{p("a"), t("1")}, {p("b"), t("2")}, {p("c"), t("1")}}, //
                     });
    }
}
