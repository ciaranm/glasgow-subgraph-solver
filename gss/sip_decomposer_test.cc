#include <gss/configuration.hh>
#include <gss/formats/csv.hh>
#include <gss/homomorphism.hh>
#include <gss/loooong.hh>
#include <gss/sip_decomposer.hh>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace gss;

using std::chrono::operator""s;
using std::make_shared;
using std::make_unique;
using std::mt19937;
using std::string;
using std::stringstream;
using std::uniform_int_distribution;
using std::uniform_real_distribution;
using std::vector;

namespace
{
    auto counting_params() -> HomomorphismParams
    {
        HomomorphismParams params;
        params.timeout = make_shared<Timeout>(0s);
        params.restarts_schedule = make_unique<NoRestartsSchedule>();
        params.count_solutions = true;
        return params;
    }

    struct Shape
    {
        bool loops, vertex_labels, edge_labels;
    };

    // A pattern made of a few small connected components, so that gluing them
    // together has something to do. Each component is a random tree with a few
    // extra edges, so it is connected.
    auto random_pattern(mt19937 & rand, const Shape & shape) -> InputGraph
    {
        uniform_int_distribution<int> number_of_components(2, 4), component_size(1, 3), coin(0, 1);
        uniform_real_distribution<double> unit(0.0, 1.0);

        vector<int> sizes;
        int total = 0;
        for (int c = number_of_components(rand); c > 0; --c) {
            sizes.push_back(component_size(rand));
            total += sizes.back();
        }

        InputGraph pattern(total, shape.vertex_labels, shape.edge_labels);
        auto edge_label = [&] { return string{shape.edge_labels ? (coin(rand) ? "r" : "g") : ""}; };

        int first = 0;
        for (auto size : sizes) {
            for (int v = first + 1; v < first + size; ++v)
                pattern.add_edge(uniform_int_distribution<int>(first, v - 1)(rand), v, edge_label());
            for (int v = first; v < first + size; ++v)
                for (int w = v + 1; w < first + size; ++w)
                    if (unit(rand) < 0.3)
                        pattern.add_edge(v, w, edge_label());
            if (shape.loops && unit(rand) < 0.2)
                pattern.add_edge(first, first, edge_label());
            first += size;
        }

        if (shape.vertex_labels)
            for (int v = 0; v < total; ++v)
                pattern.set_vertex_label(v, coin(rand) ? "a" : "b");

        return pattern;
    }

    auto random_target(mt19937 & rand, const Shape & shape) -> InputGraph
    {
        uniform_int_distribution<int> size(5, 8), coin(0, 1);
        uniform_real_distribution<double> unit(0.0, 1.0);

        int n = size(rand);
        InputGraph target(n, shape.vertex_labels, shape.edge_labels);
        for (int v = 0; v < n; ++v)
            for (int w = v; w < n; ++w)
                if ((v == w) ? (shape.loops && unit(rand) < 0.3) : unit(rand) < 0.5)
                    target.add_edge(v, w, shape.edge_labels ? (coin(rand) ? "r" : "g") : "");

        if (shape.vertex_labels)
            for (int v = 0; v < n; ++v)
                target.set_vertex_label(v, coin(rand) ? "a" : "b");

        return target;
    }
}

TEST_CASE("decomposition counts an edge plus an isolated vertex")
{
    auto pattern = read_csv(stringstream{// clang-format off
R"(a,b
c,
)"}, "pattern"); // clang-format on

    auto target = read_csv(stringstream{// clang-format off
R"(1,2
2,3
)"}, "target"); // clang-format on

    auto params = counting_params();
    auto result = solve_sip_by_decomposition(pattern, target, params);
    CHECK(result.complete);
    CHECK(result.solution_count == 4);
}

TEST_CASE("decomposition counts components that must avoid each other")
{
    // Two disjoint edges into a path on four vertices: 6 * 6 ways to embed the
    // edges independently, but only {1,2} and {3,4} are disjoint, so 2 * 2 * 2.
    auto pattern = read_csv(stringstream{// clang-format off
R"(a,b
c,d
)"}, "pattern"); // clang-format on

    auto target = read_csv(stringstream{// clang-format off
R"(1,2
2,3
3,4
)"}, "target"); // clang-format on

    auto params = counting_params();
    auto result = solve_sip_by_decomposition(pattern, target, params);
    CHECK(result.complete);
    CHECK(result.solution_count == 8);
}

TEST_CASE("decomposition agrees with the solver on random disconnected patterns")
{
    auto shape = GENERATE(Shape{false, false, false}, Shape{true, false, false}, Shape{false, true, false},
        Shape{false, false, true}, Shape{true, true, true});
    auto seed = GENERATE(range(0, 40));

    mt19937 rand(seed * 5 + shape.loops * 4 + shape.vertex_labels * 2 + shape.edge_labels);
    auto pattern = random_pattern(rand, shape);
    auto target = random_target(rand, shape);

    auto params = counting_params();
    auto expected = solve_homomorphism_problem(pattern, target, params);
    REQUIRE(expected.complete);

    auto decomposed = solve_sip_by_decomposition(pattern, target, params);
    INFO("seed " << seed << " loops " << shape.loops << " vertex labels " << shape.vertex_labels << " edge labels " << shape.edge_labels);
    CHECK(decomposed.complete);
    CHECK(decomposed.solution_count == expected.solution_count);
}

TEST_CASE("decomposition leaves non-injective counts alone")
{
    // The isolated-vertex fix-up assumes injectivity, so it must not fire here.
    auto pattern = read_csv(stringstream{// clang-format off
R"(a,b
c,
)"}, "pattern"); // clang-format on

    auto target = read_csv(stringstream{// clang-format off
R"(1,2
2,3
)"}, "target"); // clang-format on

    auto params = counting_params();
    params.injectivity = Injectivity::NonInjective;
    auto expected = solve_homomorphism_problem(pattern, target, params);
    auto decomposed = solve_sip_by_decomposition(pattern, target, params);
    CHECK(expected.solution_count == 4 * 3);
    CHECK(decomposed.solution_count == expected.solution_count);
}

TEST_CASE("decomposition refuses proof logging")
{
    auto pattern = read_csv(stringstream{"a,b\nc,\n"}, "pattern");
    auto target = read_csv(stringstream{"1,2\n2,3\n"}, "target");

    auto params = counting_params();
    params.proof_options = ProofOptions{};
    CHECK_THROWS_AS(solve_sip_by_decomposition(pattern, target, params), UnsupportedConfiguration);
}
