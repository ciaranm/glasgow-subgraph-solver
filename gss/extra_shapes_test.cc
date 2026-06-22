#include <gss/formats/csv.hh>
#include <gss/formats/input_graph.hh>
#include <gss/homomorphism.hh>

#include <catch2/catch_test_macros.hpp>

#include <memory>
#include <sstream>
#include <utility>

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

    // An extra "shape" supplemental graph: a single edge whose endpoints are
    // labelled "from" and "to". A supplemental edge between two graph vertices is
    // added when the shape embeds between them, so this one tracks plain adjacency.
    auto from_to_edge_shape() -> InputGraph
    {
        InputGraph shape{2, true, false};
        shape.add_edge(0, 1);
        shape.set_vertex_label(0, "from");
        shape.set_vertex_label(1, "to");
        return shape;
    }
}

// Extra shapes are sound supplemental filters derived from the same instance, so
// they prune the search but never change which mappings exist: the solution count
// must be identical with and without them. (Before the CLI binding was fixed, the
// --shape option was silently ignored; this guards the underlying feature.)
TEST_CASE("an extra shape preserves the solution count")
{
    auto pattern = read_csv(stringstream{"a,b\nb,c\n"}, "p");     // path a-b-c
    auto target = read_csv(stringstream{"1,2\n2,3\n3,4\n"}, "t"); // path 1-2-3-4

    auto baseline = solve_homomorphism_problem(pattern, target, make_params()).solution_count;
    CHECK(baseline > loooong{0}); // the instance is satisfiable, so the check below is meaningful

    SECTION("injective, count 1")
    {
        auto params = make_params();
        params.extra_shapes.emplace_back(make_unique<InputGraph>(from_to_edge_shape()), true, 1);
        CHECK(solve_homomorphism_problem(pattern, target, params).solution_count == baseline);
    }

    SECTION("non-injective, count 2 exercises the other code paths")
    {
        auto params = make_params();
        params.extra_shapes.emplace_back(make_unique<InputGraph>(from_to_edge_shape()), false, 2);
        CHECK(solve_homomorphism_problem(pattern, target, params).solution_count == baseline);
    }
}
