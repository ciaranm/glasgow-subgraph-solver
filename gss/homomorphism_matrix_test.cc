#include <gss/formats/csv.hh>
#include <gss/homomorphism.hh>

#include <catch2/catch_test_macros.hpp>

#include <memory>
#include <sstream>

using namespace gss;

using std::make_shared;
using std::make_unique;
using std::stringstream;

using std::chrono::operator""s;

namespace
{
    auto make_params() -> HomomorphismParams
    {
        HomomorphismParams params;
        params.timeout = make_shared<Timeout>(0s);
        params.restarts_schedule = make_unique<NoRestartsSchedule>();
        params.count_solutions = true;
        return params;
    }

    auto csv(const char * text) -> InputGraph
    {
        return read_csv(stringstream{text}, "g");
    }
}

// ---------------------------------------------------------------------------
// Induced vs non-induced
// ---------------------------------------------------------------------------

TEST_CASE("induced mapping must preserve non-edges")
{
    auto p3 = csv("a,b\nb,c\n");      // path a-b-c (a-c is a non-edge)
    auto k3 = csv("1,2\n2,3\n1,3\n"); // triangle (no non-edges)

    SECTION("non-induced: all six injective placements work")
    {
        auto params = make_params();
        CHECK(solve_homomorphism_problem(p3, k3, params).solution_count == 6);
    }

    SECTION("induced: the pattern non-edge cannot be realised in a triangle")
    {
        auto params = make_params();
        params.induced = true;
        CHECK(solve_homomorphism_problem(p3, k3, params).solution_count == 0);
    }
}

TEST_CASE("induced mapping into a graph that has non-edges")
{
    auto two_isolated = csv("a,\nb,\n"); // two vertices, no edge
    auto p3 = csv("1,2\n2,3\n");         // path: only non-adjacent pair is {1,3}

    SECTION("non-induced: any two distinct targets")
    {
        auto params = make_params();
        CHECK(solve_homomorphism_problem(two_isolated, p3, params).solution_count == 6);
    }

    SECTION("induced: only the non-adjacent target pair")
    {
        auto params = make_params();
        params.induced = true;
        CHECK(solve_homomorphism_problem(two_isolated, p3, params).solution_count == 2);
    }
}

// ---------------------------------------------------------------------------
// Injectivity modes
// ---------------------------------------------------------------------------

TEST_CASE("injectivity modes on a path-into-triangle")
{
    auto p4 = csv("a,b\nb,c\nc,d\n"); // path a-b-c-d
    auto k3 = csv("1,2\n2,3\n1,3\n"); // triangle

    SECTION("injective is impossible (four vertices into three)")
    {
        auto params = make_params();
        params.injectivity = Injectivity::Injective;
        CHECK(solve_homomorphism_problem(p4, k3, params).solution_count == 0);
    }

    SECTION("locally injective: neighbours stay distinct")
    {
        auto params = make_params();
        params.injectivity = Injectivity::LocallyInjective;
        CHECK(solve_homomorphism_problem(p4, k3, params).solution_count == 6);
    }

    SECTION("non-injective: every length-three walk")
    {
        auto params = make_params();
        params.injectivity = Injectivity::NonInjective;
        CHECK(solve_homomorphism_problem(p4, k3, params).solution_count == 24);
    }
}

// ---------------------------------------------------------------------------
// Directed
// ---------------------------------------------------------------------------

TEST_CASE("directed edges must be mapped respecting orientation")
{
    auto arc = csv("a>b\n");         // single directed edge a -> b
    auto dipath = csv("1>2\n2>3\n"); // directed path 1 -> 2 -> 3

    auto params = make_params();
    // a->b can land on 1->2 or 2->3, but not against an orientation.
    CHECK(solve_homomorphism_problem(arc, dipath, params).solution_count == 2);
}

// ---------------------------------------------------------------------------
// Labels
// ---------------------------------------------------------------------------

TEST_CASE("vertex labels constrain the mapping")
{
    auto one_red = csv("a,,red\n"); // single vertex labelled red
    auto red_and_blue = csv("1,,red\n2,,blue\n");

    auto params = make_params();
    // The lone pattern vertex can only sit on the red target.
    CHECK(solve_homomorphism_problem(one_red, red_and_blue, params).solution_count == 1);
}

TEST_CASE("edge labels constrain the mapping")
{
    auto labelled_edge = csv("a,b,x\n");     // edge a-b with label x
    auto two_labels = csv("1,2,x\n2,3,y\n"); // edge 1-2 labelled x, 2-3 labelled y

    auto params = make_params();
    // a-b can only map onto the x-labelled edge, in either direction.
    CHECK(solve_homomorphism_problem(labelled_edge, two_labels, params).solution_count == 2);
}

// ---------------------------------------------------------------------------
// Search-configuration invariants
// ---------------------------------------------------------------------------

TEST_CASE("value-ordering heuristics do not change the solution count")
{
    auto p3 = csv("a,b\nb,c\n");
    auto k3 = csv("1,2\n2,3\n1,3\n");

    for (auto heuristic : {ValueOrdering::None, ValueOrdering::Biased, ValueOrdering::Degree,
             ValueOrdering::AntiDegree, ValueOrdering::Random}) {
        auto params = make_params();
        params.value_ordering_heuristic = heuristic;
        CHECK(solve_homomorphism_problem(p3, k3, params).solution_count == 6);
    }
}
