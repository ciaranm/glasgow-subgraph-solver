#include <gss/formats/csv.hh>
#include <gss/formats/input_graph.hh>
#include <gss/homomorphism.hh>

#include <gss/configuration.hh>

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace gss;

using std::any_of;
using std::find;
using std::make_shared;
using std::make_unique;
using std::string;
using std::stringstream;
using std::vector;

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

    // A shape that does not look the same from both ends: an edge from "from" to "to", and
    // a pendant vertex on "from". It relates v and w when they are adjacent and v has some
    // other neighbour, so which of the two orientations gets tested matters. The pendant is
    // labelled "", which is what the master graph expects of an interior vertex, and which
    // is how an API caller could always write it.
    auto lopsided_shape() -> InputGraph
    {
        InputGraph shape{3, true, false};
        shape.add_edge(0, 1);
        shape.add_edge(0, 2);
        shape.set_vertex_label(0, "from");
        shape.set_vertex_label(1, "to");
        shape.set_vertex_label(2, "");
        return shape;
    }

    auto with_shape(HomomorphismParams params, InputGraph shape, bool injective = true, int count = 1) -> HomomorphismParams
    {
        params.extra_shapes.emplace_back(make_unique<InputGraph>(std::move(shape)), injective, count);
        return params;
    }

    auto stat(const HomomorphismResult & result, const string & key) -> string
    {
        for (auto & line : result.extra_stats)
            if (line.starts_with(key + " = "))
                return line.substr(key.size() + 3);
        return "";
    }

    auto words(const string & s) -> vector<string>
    {
        vector<string> result;
        stringstream ss{s};
        for (string w; ss >> w;)
            result.push_back(w);
        return result;
    }

    auto numbers_after(const string & line, const string & key) -> vector<unsigned long long>
    {
        vector<unsigned long long> result;
        for (auto & w : words(line))
            if (w.starts_with(key + ":")) {
                stringstream ss{w.substr(key.size() + 1)};
                for (string n; getline(ss, n, ',');)
                    result.push_back(stoull(n));
            }
        return result;
    }

    // Did the shape graph remove any value, either from an initial domain by its degree bound
    // or during search? Needs record_filter_activations. The initial counts are indexed by
    // graph with the original graph first, the search counts by supplemental graph only.
    auto shape_graph_removed_something(const HomomorphismResult & result) -> bool
    {
        auto names = words(stat(result, "supplemental_graph_names"));
        auto at = find(names.begin(), names.end(), "extra_shape");
        if (at == names.end())
            return false;
        auto g = unsigned(at - names.begin());
        auto initial = numbers_after(stat(result, "filter_activations_initial"), "degree");
        auto search = numbers_after(stat(result, "filter_activations_search"), "supplemental");
        return (g < initial.size() && initial[g] > 0) || (g - 1 < search.size() && search[g - 1] > 0);
    }

    auto built_a_shape_graph(const HomomorphismResult & result) -> bool
    {
        return any_of(result.extra_stats.begin(), result.extra_stats.end(), [](const string & line) {
            return line.starts_with("supplemental_graph_names =") && line.find("extra_shape") != string::npos;
        });
    }
}

// Extra shapes are sound supplemental filters derived from the same instance, so
// they prune the search but never change which mappings exist: the solution count
// must be identical with and without them. (Before the CLI binding was fixed, the
// --shape option was silently ignored; this guards the underlying feature.)
TEST_CASE("an extra shape preserves the solution count")
{
    auto pattern = read_csv(stringstream{"a,b\nb,c\n"}, "p"); // path a-b-c
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

// Issue #99: the master graph was built from the arc v -> w for w < v only, which on a
// digraph is an orientation picked by vertex numbering. Two in-stars have two mappings
// between them, and the shape used to delete one.
TEST_CASE("an extra shape preserves the solution count on directed graphs")
{
    auto pattern = read_csv(stringstream{"a>c\nb>c\n"}, "p");
    auto target = read_csv(stringstream{"1>3\n2>3\n"}, "t");

    auto params = with_shape(make_params(), from_to_edge_shape());
    auto result = solve_homomorphism_problem(pattern, target, params);
    CHECK(built_a_shape_graph(result));
    CHECK(result.solution_count == loooong{2});
}

// Issue #99: each pair was tested with "from" on the higher-numbered vertex only and then
// marked both ways, so a shape that is not symmetric in its ends deleted solutions even on
// undirected graphs. Here the path's reversal was lost.
TEST_CASE("an extra shape that is not symmetric in its ends preserves the solution count")
{
    auto pattern = read_csv(stringstream{"a,b\nb,c\n"}, "p");
    auto target = read_csv(stringstream{"1,2\n2,3\n"}, "t");

    auto params = with_shape(make_params(), lopsided_shape());
    auto result = solve_homomorphism_problem(pattern, target, params);
    CHECK(built_a_shape_graph(result));
    CHECK(result.solution_count == loooong{2});
}

// Issue #99: nothing guarded --shape, and its relation is only preserved when the mapping
// keeps the shape's vertices apart. Mapping an edge onto a single loop is a homomorphism,
// and a locally injective one, and the shape used to prune it.
TEST_CASE("an extra shape is not used without full injectivity")
{
    auto pattern = read_csv(stringstream{"a,b\n"}, "p");
    auto target = read_csv(stringstream{"1,1\n"}, "t");

    for (auto injectivity : {Injectivity::NonInjective, Injectivity::LocallyInjective}) {
        auto params = with_shape(make_params(), from_to_edge_shape());
        params.injectivity = injectivity;
        auto result = solve_homomorphism_problem(pattern, target, params);
        CHECK(! built_a_shape_graph(result));
        CHECK(result.solution_count == loooong{1});
    }
}

TEST_CASE("no supplementals means no extra shapes either")
{
    auto pattern = read_csv(stringstream{"a,b\nb,c\n"}, "p");
    auto target = read_csv(stringstream{"1,2\n2,3\n3,4\n"}, "t");

    auto params = with_shape(make_params(), from_to_edge_shape());
    params.no_supplementals = true;
    CHECK(! built_a_shape_graph(solve_homomorphism_problem(pattern, target, params)));
}

// Issue #99: every master-graph vertex other than the pair being tested is labelled "", and
// since #87 a graph file cannot give a vertex the empty label when others are labelled, so a
// shape with an interior vertex could not be read at all. Any other label now means interior.
TEST_CASE("an extra shape with an interior vertex can be read from a file")
{
    // A triangle through the two ends: it relates the pairs that are an edge of a triangle.
    auto shape = read_csv(stringstream{"f,t\nt,m\nm,f\nf,,from\nt,,to\nm,,inner\n"}, "s");

    // Every pattern vertex is on a triangle and 4 is not, although 4 has the degree for it
    // and is within distance two of everything. So with the exact-path graphs off, only the
    // shape graph can rule out mapping a pattern vertex to 4, and it should.
    auto pattern = read_csv(stringstream{"a,b\nb,c\nc,a\n"}, "p");
    auto target = read_csv(stringstream{"1,2\n2,3\n3,1\n3,4\n4,5\n"}, "t");

    auto params = with_shape(make_params(), std::move(shape));
    params.number_of_exact_path_graphs = 0;
    params.no_nds = true;
    params.record_filter_activations = true;
    auto result = solve_homomorphism_problem(pattern, target, params);
    CHECK(shape_graph_removed_something(result));
    CHECK(result.solution_count == loooong{6});
}

TEST_CASE("an extra shape must say which of its vertices are the ends")
{
    auto pattern = read_csv(stringstream{"a,b\n"}, "p");
    auto target = read_csv(stringstream{"1,2\n"}, "t");

    SECTION("unlabelled")
    {
        auto params = with_shape(make_params(), read_csv(stringstream{"f,t\n"}, "s"));
        CHECK_THROWS_AS(solve_homomorphism_problem(pattern, target, params), UnsupportedConfiguration);
    }

    SECTION("two froms")
    {
        auto params = with_shape(make_params(), read_csv(stringstream{"f,t\nf,,from\nt,,from\n"}, "s"));
        CHECK_THROWS_AS(solve_homomorphism_problem(pattern, target, params), UnsupportedConfiguration);
    }

    SECTION("edge labels")
    {
        auto params = with_shape(make_params(), read_csv(stringstream{"f,t,red\nf,,from\nt,,to\n"}, "s"));
        CHECK_THROWS_AS(solve_homomorphism_problem(pattern, target, params), UnsupportedConfiguration);
    }
}
