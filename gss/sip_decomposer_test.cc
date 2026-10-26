#include <gss/configuration.hh>
#include <gss/counting_statement.hh>
#include <gss/formats/csv.hh>
#include <gss/homomorphism.hh>
#include <gss/loooong.hh>
#include <gss/sip_decomposer.hh>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include <algorithm>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace gss;

using std::chrono::operator""s;
using std::make_shared;
using std::make_unique;
using std::map;
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

    // Every assignment of target vertices to the statement's variables, by brute force,
    // with each node's value under it. Only for tiny instances.
    template <typename F_>
    auto for_each_assignment(const CountingStatement & s, const F_ & f) -> void
    {
        int n = s.domains.at(0).size();
        vector<int> assignment(s.variables.size(), 0);
        while (true) {
            vector<bool> value(s.nodes.size());
            for (unsigned i = 0; i < s.nodes.size(); ++i) {
                auto & node = s.nodes[i];
                switch (node.kind) {
                case CountingNodeKind::Cube:
                    value[i] = true;
                    for (auto & [v, t] : node.literals)
                        value[i] = value[i] && assignment[v] == t;
                    break;
                case CountingNodeKind::False: value[i] = false; break;
                case CountingNodeKind::Product:
                    value[i] = true;
                    for (auto & c : node.children)
                        value[i] = value[i] && value[c];
                    break;
                case CountingNodeKind::Sum:
                    value[i] = false;
                    for (auto & c : node.children)
                        value[i] = value[i] || value[c];
                    break;
                case CountingNodeKind::Difference: value[i] = value[node.children[0]] && ! value[node.children[1]]; break;
                }
            }
            f(assignment, value);

            unsigned v = 0;
            for (; v < assignment.size(); ++v)
                if (++assignment[v] < n)
                    break;
                else
                    assignment[v] = 0;
            if (v == assignment.size())
                return;
        }
    }

    // A PB constraint, normalised to sum of coefficient * literal >= rhs with positive
    // coefficients, where a literal is a variable number and whether it is negated.
    struct Constraint
    {
        vector<std::tuple<long, int, bool>> terms;
        long rhs = 0;
    };

    auto parse_opb(const string & text, map<string, int> & variables) -> vector<Constraint>
    {
        vector<Constraint> result;
        stringstream in{text};
        string line;
        while (getline(in, line)) {
            if (line.empty() || line[0] == '*')
                continue;
            stringstream tokens{line};
            Constraint c;
            string word;
            long rhs = 0;
            while (tokens >> word) {
                if (word[0] == '@')
                    continue;
                if (word == ">=") {
                    tokens >> rhs;
                    break;
                }
                long coeff = std::stol(word);
                string literal;
                tokens >> literal;
                bool negated = literal[0] == '~';
                if (negated)
                    literal = literal.substr(1);
                auto id = variables.emplace(literal, variables.size()).first->second;
                if (coeff < 0) {
                    // c l = c - c ~l, so move c over to the right
                    c.rhs -= coeff;
                    coeff = -coeff;
                    negated = ! negated;
                }
                c.terms.emplace_back(coeff, id, negated);
            }
            c.rhs += rhs;
            result.push_back(c);
        }
        return result;
    }

    // Unit propagates the constraints from the given partial assignment (-1 for unset),
    // returning false on a conflict.
    auto propagate(const vector<Constraint> & constraints, vector<int> & values) -> bool
    {
        bool changed = true;
        while (changed) {
            changed = false;
            for (auto & c : constraints) {
                long slack = -c.rhs;
                for (auto & [coeff, v, negated] : c.terms)
                    if (values[v] == -1 || values[v] == (negated ? 0 : 1))
                        slack += coeff;
                if (slack < 0)
                    return false;
                for (auto & [coeff, v, negated] : c.terms)
                    if (values[v] == -1 && coeff > slack) {
                        values[v] = negated ? 0 : 1;
                        changed = true;
                    }
            }
        }
        return true;
    }

    // Checks everything about a statement that can be checked by brute force on a tiny
    // instance: its side conditions hold, its count is the number of assignments under
    // which the root holds, and the PB translation has exactly those solutions, with
    // every gate determined by propagation from the x variables.
    auto check_statement_by_brute_force(const CountingStatement & s) -> void
    {
        auto count = count_counting_statement(s);

        long root_count = 0;
        bool side_conditions_hold = true;
        for_each_assignment(s, [&](const vector<int> &, const vector<bool> & value) {
            if (value[s.root])
                ++root_count;
            for (auto & node : s.nodes)
                if (node.kind == CountingNodeKind::Sum) {
                    int holding = 0;
                    for (auto & c : node.children)
                        holding += value[c];
                    side_conditions_hold = side_conditions_hold && holding <= 1;
                }
                else if (node.kind == CountingNodeKind::Difference)
                    side_conditions_hold = side_conditions_hold && (value[node.children[0]] || ! value[node.children[1]]);
        });
        CHECK(side_conditions_hold);
        CHECK(count.count == root_count);

        stringstream opb;
        write_counting_statement_opb(opb, s);
        map<string, int> names;
        auto constraints = parse_opb(opb.str(), names);

        vector<vector<int>> x(s.variables.size());
        for (unsigned v = 0; v < s.variables.size(); ++v)
            for (auto & t : s.domains.at(0))
                x[v].push_back(names.emplace("x" + s.variables[v].first + "_" + t, names.size()).first->second);

        long opb_count = 0;
        bool all_determined = true;
        for_each_assignment(s, [&](const vector<int> & assignment, const vector<bool> &) {
            vector<int> values(names.size(), -1);
            for (unsigned v = 0; v < assignment.size(); ++v)
                for (unsigned t = 0; t < x[v].size(); ++t)
                    values[x[v][t]] = (int(t) == assignment[v]);
            if (propagate(constraints, values)) {
                for (auto & v : values)
                    all_determined = all_determined && v != -1;
                ++opb_count;
            }
        });
        CHECK(all_determined);
        CHECK(opb_count == root_count);

        // and nothing that gives a variable no value, or two, is a solution: every other
        // variable takes its first value, and v takes the values given
        bool exactly_one_enforced = true;
        auto refuted = [&](unsigned v, const vector<unsigned> & v_values) {
            vector<int> values(names.size(), -1);
            for (unsigned w = 0; w < x.size(); ++w)
                for (unsigned t = 0; t < x[w].size(); ++t)
                    values[x[w][t]] = (w == v) ? (std::find(v_values.begin(), v_values.end(), t) != v_values.end()) : (t == 0);
            return ! propagate(constraints, values);
        };
        for (unsigned v = 0; v < x.size(); ++v) {
            exactly_one_enforced = exactly_one_enforced && refuted(v, {});
            for (unsigned t = 1; t < x[v].size(); ++t)
                exactly_one_enforced = exactly_one_enforced && refuted(v, {0, t});
        }
        CHECK(exactly_one_enforced);
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

TEST_CASE("decomposition statements count what the solver counts")
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

    CountingStatement statement;
    auto decomposed = solve_sip_by_decomposition(pattern, target, params, &statement);
    INFO("seed " << seed << " loops " << shape.loops << " vertex labels " << shape.vertex_labels << " edge labels " << shape.edge_labels);
    REQUIRE(decomposed.complete);
    CHECK(decomposed.solution_count == expected.solution_count);

    // and so does what an independent counter reads back
    stringstream text;
    write_counting_statement(text, statement);
    CHECK(count_counting_statement(read_counting_statement(text)).count == expected.solution_count);
}

TEST_CASE("decomposition statements hold up under brute force")
{
    auto shape = GENERATE(Shape{false, false, false}, Shape{true, false, false}, Shape{false, true, false},
        Shape{false, false, true}, Shape{true, true, true});
    auto seed = GENERATE(range(0, 60));

    mt19937 rand(seed * 5 + shape.loops * 4 + shape.vertex_labels * 2 + shape.edge_labels);
    auto pattern = random_pattern(rand, shape);
    auto target = random_target(rand, shape);

    // brute force is over every assignment, so keep to the small ones
    long assignments = 1;
    for (int v = 0; v < pattern.size(); ++v)
        assignments *= target.size();
    if (assignments > 10000)
        return;

    auto params = counting_params();
    CountingStatement statement;
    auto decomposed = solve_sip_by_decomposition(pattern, target, params, &statement);
    REQUIRE(decomposed.complete);

    INFO("seed " << seed << " loops " << shape.loops << " vertex labels " << shape.vertex_labels << " edge labels " << shape.edge_labels);
    check_statement_by_brute_force(statement);
}

TEST_CASE("the two-edges-into-P4 statement has the expected shape")
{
    auto pattern = read_csv(stringstream{"a,b\nc,d\n"}, "pattern");
    auto target = read_csv(stringstream{"1,2\n2,3\n3,4\n"}, "target");

    auto params = counting_params();
    CountingStatement statement;
    auto result = solve_sip_by_decomposition(pattern, target, params, &statement);
    CHECK(result.solution_count == 8);

    // the root is Prod minus the sum of the six ways to glue, 36 - (4 * 4 + 2 * 6)
    auto & root = statement.nodes.at(statement.root);
    REQUIRE(root.kind == CountingNodeKind::Difference);
    auto & overlaps = statement.nodes.at(root.children[1]);
    REQUIRE(overlaps.kind == CountingNodeKind::Sum);
    CHECK(overlaps.children.size() == 6);
    check_statement_by_brute_force(statement);
}

TEST_CASE("decomposition statements need injective non-induced counting")
{
    auto pattern = read_csv(stringstream{"a,b\nc,\n"}, "pattern");
    auto target = read_csv(stringstream{"1,2\n2,3\n"}, "target");

    auto params = counting_params();
    params.induced = true;
    CountingStatement statement;
    CHECK_THROWS_AS(solve_sip_by_decomposition(pattern, target, params, &statement), UnsupportedConfiguration);
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
