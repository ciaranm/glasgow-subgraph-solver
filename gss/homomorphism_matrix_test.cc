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
    auto p3 = csv("a,b\nb,c\n"); // path a-b-c (a-c is a non-edge)
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
    auto p3 = csv("1,2\n2,3\n"); // path: only non-adjacent pair is {1,3}

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
    auto arc = csv("a>b\n"); // single directed edge a -> b
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
    auto labelled_edge = csv("a,b,x\n"); // edge a-b with label x
    auto two_labels = csv("1,2,x\n2,3,y\n"); // edge 1-2 labelled x, 2-3 labelled y

    auto params = make_params();
    // a-b can only map onto the x-labelled edge, in either direction.
    CHECK(solve_homomorphism_problem(labelled_edge, two_labels, params).solution_count == 2);
}

// A self-loop is an edge, so its label has to match like any other edge's. This is the one
// edge label the searcher's label check cannot see: loops are stripped out of the adjacency
// rows, and forward checking only ever compares a pair of distinct pattern vertices, so it
// is the loop-compatibility check that has to do it (issue #92).
TEST_CASE("edge labels on self-loops constrain the mapping")
{
    auto red_loop = csv("a,a,red\n");

    auto params = make_params();
    CHECK(solve_homomorphism_problem(red_loop, csv("1,1,blue\n"), params).solution_count == 0);
    CHECK(solve_homomorphism_problem(red_loop, csv("1,1,red\n"), params).solution_count == 1);
}

// ---------------------------------------------------------------------------
// Clique-size constraints
// ---------------------------------------------------------------------------

// The clique-size filter (--cliques) maps a pattern vertex in a k-clique only to a target
// vertex in a k-clique, which needs the k pattern vertices to reach k distinct targets.
// Injectivity gives that, and so -- on the original graph pair only -- does a loopless
// target, since collapsing two adjacent vertices would need a self-loop on the image. Both
// arguments had been taken for granted rather than checked, so the filter silently deleted
// solutions in the two cases below (issue #91).
TEST_CASE("clique-size constraints do not change the solution count")
{
    SECTION("non-injective into a target with a loop")
    {
        auto edge = csv("a,b\n");
        auto oneloop = csv("1,1\n"); // one vertex, self-loop

        auto params = make_params();
        params.injectivity = Injectivity::NonInjective;
        params.clique_size_constraints = true;
        // both pattern vertices onto the loop
        CHECK(solve_homomorphism_problem(edge, oneloop, params).solution_count == 1);
    }

    SECTION("non-injective with clique sizes on the supplemental graphs, no loops anywhere")
    {
        auto star = csv("b,a\nb,c\nb,d\n"); // K_{1,3}, centre b
        auto edge = csv("1,2\n");

        auto params = make_params();
        params.injectivity = Injectivity::NonInjective;
        params.clique_size_constraints = true;
        params.clique_size_constraints_on_supplementals = true;
        // the three leaves are a triangle in the pattern's distance-2 graph, but they may
        // legitimately share an image: centre on either end of the edge, leaves on the other
        CHECK(solve_homomorphism_problem(star, edge, params).solution_count == 2);
    }

    SECTION("injective: the filter still fires, and still counts the same")
    {
        auto triangle = csv("a,b\nb,c\na,c\n");
        auto triangle_and_c4 = csv("1,2\n2,3\n1,3\n4,5\n5,6\n6,7\n7,4\n");

        // supplementals and NDS off, so that what is left to prune the triangle-free half of
        // the target is the clique-size filter itself
        auto without = make_params();
        without.no_supplementals = true;
        without.no_nds = true;
        auto without_result = solve_homomorphism_problem(triangle, triangle_and_c4, without);

        auto with = make_params();
        with.no_supplementals = true;
        with.no_nds = true;
        with.clique_size_constraints = true;
        auto with_result = solve_homomorphism_problem(triangle, triangle_and_c4, with);

        CHECK(with_result.solution_count == without_result.solution_count);
        CHECK(with_result.solution_count == 6);
        CHECK(with_result.nodes < without_result.nodes);
    }
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
