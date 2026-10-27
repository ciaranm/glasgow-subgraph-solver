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
#include <filesystem>
#include <map>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <unistd.h>

using namespace gss;

using std::chrono::operator""s;
using std::getline;
using std::make_shared;
using std::make_unique;
using std::map;
using std::mt19937;
using std::optional;
using std::pair;
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

    // A number of hubs to force, or -1 to let the decomposer choose.
    auto hubs_option(int hubs) -> optional<int>
    {
        return hubs == -1 ? std::nullopt : optional<int>{hubs};
    }

    auto has_stat(const HomomorphismResult & result, const string & stat) -> bool
    {
        return std::find(result.extra_stats.begin(), result.extra_stats.end(), stat) != result.extra_stats.end();
    }

    // Writes a statement, and reads it back, which checks that it is well-formed.
    auto statement_for(const InputGraph & pattern, const InputGraph & target, const HomomorphismParams & params,
        HomomorphismResult & result, optional<int> hubs = std::nullopt, bool check_whole_instance = true) -> optional<CountingStatement>
    {
        stringstream out;
        result = solve_sip_by_decomposition(pattern, target, params, DecompositionOptions{&out, hubs, check_whole_instance});
        if (out.str().empty())
            return std::nullopt;
        return read_counting_statement(out);
    }

    // What each node says about an assignment f, straight from the definitions in
    // dev_docs/counting-statements.md.
    auto evaluate(const CountingStatement & statement, const vector<int> & f) -> vector<bool>
    {
        vector<bool> result;
        for (auto & node : statement.nodes) {
            bool value = false;
            switch (node.kind) {
            case CountingNodeKind::Cube:
                value = true;
                for (auto & [u, t] : node.assignment)
                    value = value && f[u] == t;
                break;
            case CountingNodeKind::And:
                value = true;
                for (auto & c : node.children)
                    value = value && result[c];
                break;
            case CountingNodeKind::Or:
                for (auto & c : node.children)
                    value = value || result[c];
                break;
            case CountingNodeKind::Minus:
                value = result[node.children[0]] && ! result[node.children[1]];
                break;
            }
            result.push_back(value);
        }
        return result;
    }

    // A row of one of the translator's OPB files, which only ever use >= or <= and
    // integer coefficients on plain or negated variables. Variables are numbered by
    // name as they are first seen, so that evaluating a row is cheap.
    struct Row
    {
        vector<pair<long, pair<bool, int>>> terms;
        bool at_most;
        long degree;

        auto satisfied(const vector<char> & values) const -> bool
        {
            long lhs = 0;
            for (auto & [coefficient, literal] : terms)
                if (bool(values[literal.second]) != literal.first)
                    lhs += coefficient;
            return at_most ? lhs <= degree : lhs >= degree;
        }

        auto mentions(int variable) const -> bool
        {
            for (auto & [_, literal] : terms)
                if (literal.second == variable)
                    return true;
            return false;
        }
    };

    auto parse_rows(const string & text, map<string, int> & names, vector<string> & preserved) -> vector<Row>
    {
        auto number = [&](const string & name) { return names.emplace(name, names.size()).first->second; };
        vector<Row> rows;
        stringstream in{text};
        string line;
        while (getline(in, line)) {
            stringstream words{line};
            string word;
            words >> word;
            if (word == "preserved:") {
                while (words >> word && word != ";")
                    preserved.push_back(word);
                continue;
            }
            Row row;
            while (true) {
                if (word == ">=" || word == "<=") {
                    row.at_most = word == "<=";
                    words >> row.degree;
                    break;
                }
                long coefficient = std::stol(word);
                string literal;
                words >> literal;
                bool negated = literal[0] == '~';
                row.terms.emplace_back(coefficient, pair{negated, number(negated ? literal.substr(1) : literal)});
                words >> word;
            }
            rows.push_back(row);
        }
        return rows;
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

    // every number of hubs gives the same count; -1 is the decomposer's own choice. The
    // checks on the whole instance rule out many of these, so the decomposition itself
    // also runs without them.
    for (bool check : {true, false})
        for (int hubs : {-1, 0, 1, 2, 3}) {
            auto decomposed = solve_sip_by_decomposition(pattern, target, params, DecompositionOptions{nullptr, hubs_option(hubs), check});
            INFO("seed " << seed << " loops " << shape.loops << " vertex labels " << shape.vertex_labels << " edge labels " << shape.edge_labels << " hubs " << hubs << " check " << check);
            CHECK(decomposed.complete);
            CHECK(decomposed.solution_count == expected.solution_count);
            if (hubs != -1)
                CHECK(has_stat(decomposed, "decomposition_hubs = " + std::to_string(hubs)));
            if (has_stat(decomposed, "decomposition_rejected_before_search = true"))
                CHECK((check && expected.solution_count == 0));
        }
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

TEST_CASE("decomposition statements count what the solver counts on random disconnected patterns")
{
    auto shape = GENERATE(Shape{false, false, false}, Shape{true, false, false}, Shape{false, true, false},
        Shape{false, false, true}, Shape{true, true, true});
    auto seed = GENERATE(range(0, 40));

    // A statement has a node for every quotient, not just one per isomorphism class,
    // so the biggest of the patterns above can take tens of seconds.
    mt19937 rand(seed * 5 + shape.loops * 4 + shape.vertex_labels * 2 + shape.edge_labels);
    auto pattern = [&] {
        while (true)
            if (auto p = random_pattern(rand, shape); p.size() <= 9)
                return p;
    }();
    auto target = random_target(rand, shape);

    auto params = counting_params();
    auto expected = solve_homomorphism_problem(pattern, target, params);
    REQUIRE(expected.complete);

    for (bool check : {true, false})
        for (int hubs : {-1, 0, 1, 2, 3}) {
            HomomorphismResult result;
            auto statement = statement_for(pattern, target, params, result, hubs_option(hubs), check);
            INFO("seed " << seed << " loops " << shape.loops << " vertex labels " << shape.vertex_labels << " edge labels " << shape.edge_labels << " hubs " << hubs << " check " << check);
            CHECK(result.complete);
            CHECK(result.solution_count == expected.solution_count);
            if (statement) {
                CHECK(statement->domain_sizes == vector<int>(pattern.size(), target.size()));
                CHECK(count_counting_statement(*statement) == expected.solution_count);
            }
            else
                CHECK(expected.solution_count == 0);
        }
}

TEST_CASE("decomposition statements mean what they say by brute force")
{
    // Every assignment of a tiny instance against the semantics of each node and the
    // side conditions, and the solutions plus a sample of the other assignments against
    // the translator's rows.
    auto shape = GENERATE(Shape{false, false, false}, Shape{true, false, false}, Shape{false, true, false},
        Shape{false, false, true}, Shape{true, true, true});
    auto seed = GENERATE(range(0, 30));
    auto hubs = GENERATE(-1, 0, 1, 2);
    // without the checks on the whole instance, which would leave many with no statement

    // no more pattern vertices than target vertices, or the pigeonhole would mean that
    // most of these write no statement at all
    mt19937 rand(1000 + seed * 5 + shape.loops * 4 + shape.vertex_labels * 2 + shape.edge_labels);
    auto pattern = [&] {
        while (true)
            if (auto p = random_pattern(rand, shape); p.size() <= 5)
                return p;
    }();
    auto target = [&] {
        while (true)
            if (auto t = random_target(rand, shape); t.size() <= 5)
                return t;
    }();

    auto params = counting_params();
    auto expected = solve_homomorphism_problem(pattern, target, params);
    REQUIRE(expected.complete);

    HomomorphismResult result;
    auto statement = statement_for(pattern, target, params, result, hubs_option(hubs), false);
    INFO("seed " << seed << " loops " << shape.loops << " vertex labels " << shape.vertex_labels << " edge labels " << shape.edge_labels << " hubs " << hubs);
    if (! statement) {
        CHECK(expected.solution_count == 0);
        return;
    }

    REQUIRE(statement->domain_sizes == vector<int>(pattern.size(), target.size()));
    int m = pattern.size(), n = target.size(), size = statement->nodes.size();
    auto side_conditions = counting_side_conditions(*statement);

    stringstream opb, side;
    write_counting_opb(opb, *statement);
    bool has_side = write_side_condition_opb(side, *statement);

    map<string, int> names;
    vector<int> x(m * n), g(size), v(size, -1);
    for (int u = 0; u < m; ++u)
        for (int t = 0; t < n; ++t)
            x[u * n + t] = names.emplace("y" + std::to_string(u) + "_" + std::to_string(t), names.size()).first->second;
    for (int i = 0; i < size; ++i)
        g[i] = names.emplace("g" + std::to_string(i), names.size()).first->second;
    for (auto & i : side_conditions)
        v[i] = names.emplace("v" + std::to_string(i), names.size()).first->second;

    vector<string> preserved, side_preserved;
    auto rows = parse_rows(opb.str(), names, preserved);
    auto side_rows = has_side ? parse_rows(side.str(), names, side_preserved) : vector<Row>{};
    REQUIRE(names.size() == x.size() + g.size() + side_conditions.size());
    CHECK(int(preserved.size()) == m * n);
    CHECK(side_preserved.empty());

    // which rows define which node: two per variable for the exactly-one groups, then
    // each node's rows in order, then the root
    vector<pair<int, int>> rows_of_node;
    int next_row = 2 * m;
    for (auto & node : statement->nodes) {
        bool constant = node.kind != CountingNodeKind::Or && node.kind != CountingNodeKind::Minus && node.assignment.empty() && node.children.empty();
        int k = node.kind == CountingNodeKind::Minus ? 3 : constant ? 1
                                                                    : 2;
        rows_of_node.emplace_back(next_row, next_row + k);
        next_row += k;
    }
    REQUIRE(next_row + 1 == int(rows.size()));

    vector<vector<int>> side_rows_of_selector(size);
    for (int r = 0; r + 1 < int(side_rows.size()); ++r)
        for (auto & i : side_conditions)
            if (side_rows[r].mentions(v[i]))
                side_rows_of_selector[i].push_back(r);

    long total = 0, disjointness_failures = 0, containment_failures = 0, row_failures = 0, root_row_failures = 0,
         undetermined_gates = 0, side_failures = 0, rows_checked_for = 0;
    std::uniform_int_distribution<int> sample(0, 99);
    vector<int> f(m, 0);
    vector<char> values(names.size());
    while (true) {
        auto value = evaluate(*statement, f);
        total += value[statement->root];

        for (int i = 0; i < size; ++i) {
            auto & node = statement->nodes[i];
            if (node.kind == CountingNodeKind::Or) {
                int true_children = 0;
                for (auto & c : node.children)
                    true_children += value[c];
                disjointness_failures += true_children > 1;
            }
            else if (node.kind == CountingNodeKind::Minus)
                containment_failures += value[node.children[1]] && ! value[node.children[0]];
        }

        if (value[statement->root] || sample(rand) < 2) {
            ++rows_checked_for;
            for (int u = 0; u < m; ++u)
                for (int t = 0; t < n; ++t)
                    values[x[u * n + t]] = f[u] == t;
            for (int i = 0; i < size; ++i)
                values[g[i]] = value[i];
            for (auto & i : side_conditions)
                values[v[i]] = false;

            // with the gates set by the semantics, every definition row holds, the root
            // row holds exactly when the root does, and flipping any one gate breaks one
            // of its own rows, so the gates are determined
            for (int r = 0; r + 1 < int(rows.size()); ++r)
                row_failures += ! rows[r].satisfied(values);
            root_row_failures += rows.back().satisfied(values) != value[statement->root];
            for (int i = 0; i < size; ++i) {
                values[g[i]] = ! values[g[i]];
                bool broken = false;
                for (int r = rows_of_node[i].first; r < rows_of_node[i].second; ++r)
                    broken = broken || ! rows[r].satisfied(values);
                undetermined_gates += ! broken;
                values[g[i]] = ! values[g[i]];
            }

            // the side-condition instance: its definitions hold, the selector row fails
            // with every selector off, and no selector can be switched on
            if (has_side) {
                for (int r = 0; r + 1 < int(side_rows.size()); ++r)
                    side_failures += ! side_rows[r].satisfied(values);
                side_failures += side_rows.back().satisfied(values);
                for (auto & i : side_conditions) {
                    values[v[i]] = true;
                    bool broken = false;
                    for (auto & r : side_rows_of_selector[i])
                        broken = broken || ! side_rows[r].satisfied(values);
                    side_failures += ! broken;
                    values[v[i]] = false;
                }
            }
        }

        int u = 0;
        while (u < m && ++f[u] == n)
            f[u++] = 0;
        if (u == m)
            break;
    }

    CHECK(disjointness_failures == 0);
    CHECK(containment_failures == 0);
    CHECK(row_failures == 0);
    CHECK(root_row_failures == 0);
    CHECK(undetermined_gates == 0);
    CHECK(side_failures == 0);
    CHECK(rows_checked_for > 0);
    CHECK(loooong{total} == expected.solution_count);
    CHECK(count_counting_statement(*statement) == expected.solution_count);
}

TEST_CASE("decomposition statements are refused where they cannot be written")
{
    auto pattern = read_csv(stringstream{"a,b\nc,\n"}, "pattern");
    auto target = read_csv(stringstream{"1,2\n2,3\n"}, "target");

    auto params = counting_params();
    params.injectivity = Injectivity::NonInjective;
    stringstream out;
    CHECK_THROWS_AS(solve_sip_by_decomposition(pattern, target, params, DecompositionOptions{&out, std::nullopt}), UnsupportedConfiguration);
}

TEST_CASE("decomposition statements are not written when a component cannot be placed")
{
    // a triangle cannot go into a path, so there is nothing to subtract from
    auto pattern = read_csv(stringstream{"a,b\nb,c\nc,a\nd,e\n"}, "pattern");
    auto target = read_csv(stringstream{"1,2\n2,3\n3,4\n"}, "target");

    auto params = counting_params();
    HomomorphismResult result;
    CHECK(! statement_for(pattern, target, params, result));
    CHECK(result.complete);
    CHECK(result.solution_count == 0);
}

TEST_CASE("decomposition statements are not written when the pattern cannot fit by pigeonhole")
{
    // two disjoint edges and an isolated vertex need five target vertices
    auto pattern = read_csv(stringstream{"a,b\nc,d\ne,\n"}, "pattern");
    auto target = read_csv(stringstream{"1,2\n2,3\n3,4\n"}, "target");

    auto params = counting_params();
    HomomorphismResult result;
    CHECK(! statement_for(pattern, target, params, result));
    CHECK(result.complete);
    CHECK(result.solution_count == 0);
    CHECK(std::find(result.extra_stats.begin(), result.extra_stats.end(),
              "decomposition_statement = none, because the pattern has more vertices of some label than the target") != result.extra_stats.end());
}

TEST_CASE("decomposition does not glue a pattern that is bigger than its target")
{
    // without the pigeonhole, this explores every way of gluing twelve edges together
    string edges;
    for (int i = 0; i < 12; ++i)
        edges += "a" + std::to_string(i) + ",b" + std::to_string(i) + "\n";
    auto pattern = read_csv(stringstream{edges}, "pattern");
    auto target = read_csv(stringstream{"1,2\n2,3\n3,4\n4,5\n5,6\n"}, "target");

    auto params = counting_params();
    auto result = solve_sip_by_decomposition(pattern, target, params);
    CHECK(result.complete);
    CHECK(result.solution_count == 0);
    CHECK(std::find(result.extra_stats.begin(), result.extra_stats.end(), "decomposition_coarsenings_considered = 0") != result.extra_stats.end());
}

TEST_CASE("decomposition starts with the checks on the whole instance")
{
    // each edge fits on its own, but the degree sequences rule out both together
    auto pattern = read_csv(stringstream{"a,b\nc,d\n"}, "pattern");
    auto target = read_csv(stringstream{"1,2\n2,3\n4,\n"}, "target");
    auto params = counting_params();

    auto checked = solve_sip_by_decomposition(pattern, target, params);
    CHECK(checked.complete);
    CHECK(checked.solution_count == 0);
    CHECK(has_stat(checked, "decomposition_rejected_before_search = true"));
    CHECK(has_stat(checked, "decomposition_coarsenings_considered = 0"));

    HomomorphismResult with_statement;
    CHECK(! statement_for(pattern, target, params, with_statement));
    CHECK(has_stat(with_statement, "decomposition_statement = none, because the checks before search rule out every embedding"));

    auto unchecked = solve_sip_by_decomposition(pattern, target, params, DecompositionOptions{nullptr, std::nullopt, false});
    CHECK(unchecked.complete);
    CHECK(unchecked.solution_count == 0);
    CHECK(! has_stat(unchecked, "decomposition_rejected_before_search = true"));
    CHECK(! has_stat(unchecked, "decomposition_coarsenings_considered = 0"));
}

TEST_CASE("decomposition statements can say zero when the components only overlap")
{
    // two disjoint edges into a star with three leaves: each edge fits, but never both,
    // and there are enough target vertices that the pigeonhole does not see it
    auto pattern = read_csv(stringstream{"a,b\nc,d\n"}, "pattern");
    auto target = read_csv(stringstream{"1,2\n1,3\n1,4\n"}, "target");

    auto params = counting_params();
    HomomorphismResult result;
    auto statement = statement_for(pattern, target, params, result);
    REQUIRE(statement);
    CHECK(count_counting_statement(*statement) == 0);
    CHECK(result.solution_count == 0);
}

TEST_CASE("decomposition writes a proof alongside a statement")
{
    // two disjoint edges into a path on four vertices, with a proof; whether the proof
    // verifies is checked by the decomposition_proof_sweep test, when VeriPB is around
    auto pattern = read_csv(stringstream{"a,b\nc,d\n"}, "pattern");
    auto target = read_csv(stringstream{"1,2\n2,3\n3,4\n"}, "target");

    auto base = std::filesystem::temp_directory_path() / ("gss_decomposition_proof_" + std::to_string(::getpid()));
    auto params = counting_params();
    params.proof_options = ProofOptions{.opb_file = base.string() + ".opb", .log_file = base.string() + ".pbp", .super_extra_verbose = false};
    HomomorphismResult result;
    auto statement = statement_for(pattern, target, params, result);
    REQUIRE(statement);
    CHECK(result.solution_count == 8);
    CHECK(count_counting_statement(*statement) == 8);
    CHECK(std::filesystem::file_size(base.string() + ".opb") > 0);
    CHECK(std::filesystem::file_size(base.string() + ".pbp") > 0);
    std::filesystem::remove(base.string() + ".opb");
    std::filesystem::remove(base.string() + ".pbp");
}

TEST_CASE("decomposition refuses proof logging on labelled graphs")
{
    auto pattern = read_csv(stringstream{"a,b,x\nc,d,y\n"}, "pattern");
    auto target = read_csv(stringstream{"1,2,x\n2,3,y\n"}, "target");

    auto params = counting_params();
    params.proof_options = ProofOptions{};
    stringstream out;
    CHECK_THROWS_AS(solve_sip_by_decomposition(pattern, target, params, DecompositionOptions{&out, std::nullopt}), UnsupportedConfiguration);
}

TEST_CASE("decomposition refuses proof logging")
{
    auto pattern = read_csv(stringstream{"a,b\nc,\n"}, "pattern");
    auto target = read_csv(stringstream{"1,2\n2,3\n"}, "target");

    auto params = counting_params();
    params.proof_options = ProofOptions{};
    CHECK_THROWS_AS(solve_sip_by_decomposition(pattern, target, params), UnsupportedConfiguration);
}

TEST_CASE("decomposition guesses a star's centre onto a hub")
{
    // A star with four leaves into a star with 199 leaves, which a single search cannot
    // count without visiting every one of its 199 * 198 * 197 * 196 solutions. With the
    // hub guessed, the leaves are four isolated vertices with lists.
    string star;
    for (int i = 1; i <= 199; ++i)
        star += "0," + std::to_string(i) + "\n";
    auto target = read_csv(stringstream{star}, "target");
    auto pattern = read_csv(stringstream{"c,a\nc,b\nc,d\nc,e\n"}, "pattern");

    auto params = counting_params();
    loooong expected = loooong{199} * 198 * 197 * 196;
    auto result = solve_sip_by_decomposition(pattern, target, params);
    CHECK(result.complete);
    CHECK(result.solution_count == expected);
    CHECK(has_stat(result, "decomposition_hubs = 1"));

    HomomorphismResult with_statement;
    auto statement = statement_for(pattern, target, params, with_statement);
    REQUIRE(statement);
    CHECK(with_statement.solution_count == expected);
    CHECK(count_counting_statement(*statement) == expected);
    CHECK(has_stat(with_statement, "decomposition_hubs = 1"));
}

TEST_CASE("decomposition hubs are refused where they cannot be used")
{
    auto pattern = read_csv(stringstream{"a,b\nc,d\n"}, "pattern");
    auto target = read_csv(stringstream{"1,2\n2,3\n3,4\n"}, "target");
    auto params = counting_params();

    // more hubs than target vertices, or fewer than none
    CHECK_THROWS_AS(solve_sip_by_decomposition(pattern, target, params, DecompositionOptions{nullptr, 5}), UnsupportedConfiguration);
    CHECK_THROWS_AS(solve_sip_by_decomposition(pattern, target, params, DecompositionOptions{nullptr, -1}), UnsupportedConfiguration);

    // not when there is nothing to count by components
    auto non_injective = counting_params();
    non_injective.injectivity = Injectivity::NonInjective;
    CHECK_THROWS_AS(solve_sip_by_decomposition(pattern, target, non_injective, DecompositionOptions{nullptr, 1}), UnsupportedConfiguration);

    // not into a directed target
    InputGraph directed(3, false, false, true);
    directed.add_directed_edge(0, 1, "");
    directed.add_directed_edge(1, 2, "");
    CHECK_THROWS_AS(solve_sip_by_decomposition(pattern, directed, params, DecompositionOptions{nullptr, 1}), UnsupportedConfiguration);
}
