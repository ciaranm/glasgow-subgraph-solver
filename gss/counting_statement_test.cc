#include <gss/counting_statement.hh>
#include <gss/counting_statement_views.hh>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

using namespace gss;

using std::map;
using std::set;
using std::string;
using std::stringstream;
using std::vector;

namespace
{
    auto read(const string & text) -> CountingStatement
    {
        stringstream in{text};
        return read_counting_statement(in);
    }

    auto statement(const string & domain, const string & nodes, const string & root) -> string
    {
        return "((version 3) (domain " + domain + ") (nodes " + nodes + ") (root " + root + "))";
    }

    // one variable, three values, either of the first two
    const string one_of_two = statement("3", "(0 1 cube (0 0)) (1 1 cube (0 1)) (2 2 or 0 1)", "2");

    // The number of solutions of an OPB file, projected onto its preserved variables, by
    // trying every assignment: only for tiny files, in the forms write_tagged_opb writes.
    auto count_solutions(const string & opb) -> long
    {
        struct Row
        {
            vector<std::tuple<int, string, bool>> literals;
            bool at_most;
            int degree;
        };
        vector<string> preserved;
        vector<Row> rows;
        set<string> variables;
        stringstream in{opb};
        string line;
        while (std::getline(in, line)) {
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
                int coefficient = std::stoi(word);
                REQUIRE(coefficient >= 1);
                words >> word;
                bool negated = word[0] == '~';
                auto name = negated ? word.substr(1) : word;
                row.literals.emplace_back(coefficient, name, negated);
                variables.insert(name);
                words >> word;
            }
            rows.push_back(row);
        }

        vector<string> names{variables.begin(), variables.end()};
        REQUIRE(names.size() <= 20);
        set<vector<bool>> projections;
        for (unsigned long bits = 0; bits < (1ul << names.size()); ++bits) {
            map<string, bool> value;
            for (unsigned i = 0; i < names.size(); ++i)
                value[names[i]] = (bits >> i) & 1;
            bool ok = true;
            for (auto & row : rows) {
                int sum = 0;
                for (auto & [coefficient, name, negated] : row.literals)
                    sum += coefficient * (value[name] != negated);
                if (row.at_most ? sum > row.degree : sum < row.degree) {
                    ok = false;
                    break;
                }
            }
            if (ok) {
                vector<bool> projection;
                for (auto & p : preserved)
                    projection.push_back(value[p]);
                projections.insert(projection);
            }
        }
        return projections.size();
    }
}

TEST_CASE("counting statements read and count")
{
    auto s = read(one_of_two);
    CHECK(s.domain_sizes == vector<int>{3});
    CHECK(s.nodes.size() == 3);
    CHECK(s.nodes[s.root].count == 2);
}

TEST_CASE("counting statements multiply and add")
{
    // x0 in {0, 1}, x1 in {0, 1, 2}
    auto s = read(statement("3 3",
        "(0 1 cube (0 0)) (1 1 cube (0 1)) (2 2 or 0 1) "
        "(3 1 cube (1 0)) (4 1 cube (1 1)) (5 1 cube (1 2)) (6 3 or 3 4 5) "
        "(7 6 and 2 6)",
        "7"));
    CHECK(s.nodes[7].count == 6);
}

TEST_CASE("counting statements count an or's children separately, even where they overlap")
{
    auto s = read(statement("2", "(0 1 cube (0 0)) (1 1 cube (0 0)) (2 1 cube (0 1)) (3 3 or 0 1 2)", "3"));
    CHECK(s.nodes[3].count == 3);
}

TEST_CASE("counting statements allow empty cubes and ands")
{
    CHECK(read(statement("", "(0 1 cube)", "0")).nodes[0].count == 1);
    CHECK(read(statement("", "(0 1 and)", "0")).nodes[0].count == 1);
}

TEST_CASE("counting statements reject anything malformed")
{
    auto bad = GENERATE(
        // lexical and grammatical
        string{""},
        string{"("},
        one_of_two + " extra",
        string{"((version 1) (domain 1) (nodes (0 1 cube (0 0))) (root 0))"},
        // version 2 had no distinct nodes, and is not read as 3
        string{"((version 2) (domain 1) (nodes (0 1 cube (0 0))) (root 0))"},
        string{"((domain 1) (version 3) (nodes (0 1 cube (0 0))) (root 0))"},
        string{"((version 3) (domain 1) (nodes (0 1 cube (0 0))))"},
        statement("03", "(0 1 cube (0 0))", "0"),
        statement("-3", "(0 1 cube (0 0))", "0"),
        statement("3", "(0 cube (0 0))", "0"),
        statement("3", "(0 01 cube (0 0))", "0"),
        statement("3", "(0 -1 cube (0 0))", "0"),
        statement("3", "(0 1 cube (0 0)) (1 1 xor 0)", "1"),
        statement("3", "(0 1 cube (0 0)) (1 1 cube (0 1)) (2 1 minus 0 1)", "2"),
        statement("3", "(0 1 cube (0 0) (0))", "0"),
        statement("3", "(0 1 cube (0 0 1))", "0"),
        statement("3", "(0 1 cube 0 0)", "0"),
        statement("3", "(0 1 or (0 0))", "0"),
        statement("99999999999", "(0 1 cube (0 0))", "0"),
        // ids and references
        statement("3", "(1 1 cube (0 0))", "1"),
        statement("3", "(0 1 cube (0 0)) (0 1 cube (0 1))", "0"),
        statement("3", "(0 1 or 0)", "0"),
        statement("3", "(0 1 or 1) (1 1 cube (0 0))", "0"),
        statement("3", "(0 1 cube (0 0))", "1"),
        // children in strictly increasing order
        statement("3", "(0 1 cube (0 0)) (1 1 cube (0 1)) (2 2 or 1 0)", "2"),
        statement("3", "(0 1 cube (0 0)) (1 2 or 0 0)", "1"),
        statement("3 3", "(0 1 cube (0 0)) (1 1 cube (1 0)) (2 1 and 1 0)", "2"),
        // counts
        statement("3", "(0 2 cube (0 0))", "0"),
        statement("3", "(0 1 cube (0 0)) (1 1 cube (0 1)) (2 1 or 0 1)", "2"),
        statement("3 3", "(0 1 cube (0 0)) (1 1 cube (0 1)) (2 2 or 0 1) (3 1 cube (1 0)) (4 3 and 2 3)", "4"),
        // ranges
        statement("3", "(0 1 cube (1 0))", "0"),
        statement("3", "(0 1 cube (0 3))", "0"),
        // cubes
        statement("3 3", "(0 1 cube (1 0) (0 0))", "0"),
        statement("3 3", "(0 1 cube (0 0) (0 1))", "0"),
        // scopes
        statement("3 3", "(0 1 cube (0 0)) (1 1 cube (0 1) (1 1)) (2 1 and 0 1)", "2"),
        statement("3 3", "(0 1 cube (0 0) (1 0)) (1 1 cube (0 1)) (2 2 or 0 1)", "2"),
        statement("3 3", "(0 1 cube (0 0) (1 0)) (1 0 or)", "0"),
        statement("", "(0 1 cube) (1 1 and) (2 2 or 0 1)", "2"),
        // the root must cover every variable
        statement("3 3", "(0 1 cube (0 0))", "0"),
        statement("3 3", "(0 1 cube (0 0) (1 0)) (1 1 cube (0 0))", "1"),
        // each variable has its own domain
        statement("2 3", "(0 1 cube (0 2) (1 0))", "0"),
        statement("3 0", "(0 1 cube (0 0))", "0"),
        // distinct nodes
        statement("3 3", "(0 6 distinct (0 1) (0 1 2))", "1"),
        statement("3 3", "(0 5 distinct (0 1) (0 1 2))", "0"),
        statement("3 3", "(0 6 distinct (1 0) (0 1 2))", "0"),
        statement("3 3", "(0 6 distinct (0 1) (0 2 1))", "0"),
        statement("3 3", "(0 2 distinct (0 1) (0 0 1))", "0"),
        statement("3 3", "(0 0 distinct (0 1) (0))", "0"),
        statement("3 3", "(0 1 distinct () (0 1))", "0"),
        statement("3 3", "(0 2 distinct (0 1) (0 1)) (1 1 cube (0 0)) (2 2 and 0 1)", "2"),
        statement("3 2", "(0 6 distinct (0 1) (0 1 2))", "0"),
        statement("3 3", "(0 6 distinct (0 2) (0 1 2))", "0"),
        statement("3 3", "(0 6 distinct (0 1) (0 1 2) (3))", "0"),
        statement("3 3", "(0 6 distinct 0 1 (0 1 2))", "0"),
        statement("3 3", "(0 1 cube (0 0)) (1 6 distinct (0 1) (0 1 2) 0)", "1"));

    INFO(bad);
    CHECK_THROWS_AS(read(bad), CountingStatementError);
}

TEST_CASE("counting statements read counts of any size")
{
    // eighty two-valued variables, each free: a count past 64 bits, read and checked
    // exactly
    string domain = "2";
    string nodes = "(0 1 cube (0 0)) (1 1 cube (0 1)) (2 2 or 0 1)";
    int last = 2;
    loooong count = 2;
    for (int u = 1; u < 80; ++u) {
        domain += " 2";
        int a = last + 1, b = last + 2, c = last + 3, d = last + 4;
        count *= 2;
        stringstream more;
        more << " (" << a << " 1 cube (" << u << " 0)) (" << b << " 1 cube (" << u << " 1)) (" << c << " 2 or " << a << " " << b
             << ") (" << d << " " << count << " and " << last << " " << c << ")";
        nodes += more.str();
        last = d;
    }
    auto s = read(statement(domain, nodes, std::to_string(last)));
    stringstream shown;
    shown << s.nodes[last].count;
    CHECK(shown.str() == "1208925819614629174706176");
}

TEST_CASE("counting statements translate exactly as specified")
{
    // node 2 overlaps itself through node 1, node 3 is never used, node 4 shares node 0
    auto s = read(statement("2 2",
        "(0 1 cube (0 0) (1 1)) (1 1 cube (0 1) (1 0)) (2 2 or 0 1) (3 1 and) (4 3 or 0 2)",
        "4"));

    stringstream opb;
    write_tagged_opb(opb, s);
    CHECK(opb.str() ==
        "preserved: y0_0 y0_1 y1_0 y1_1 t2_0 t2_1 t4_0 t4_1 ;\n"
        "1 y0_0 1 y0_1 >= 1 ;\n"
        "1 y0_0 1 y0_1 <= 1 ;\n"
        "1 y1_0 1 y1_1 >= 1 ;\n"
        "1 y1_0 1 y1_1 <= 1 ;\n"
        "1 a4 >= 1 ;\n"
        "1 ~a0 1 t2_0 1 t4_0 >= 1 ;\n"
        "1 a0 1 ~t2_0 >= 1 ;\n"
        "1 a0 1 ~t4_0 >= 1 ;\n"
        "1 ~a0 1 y0_0 >= 1 ;\n"
        "1 ~a0 1 y1_1 >= 1 ;\n"
        "1 ~a1 1 t2_1 >= 1 ;\n"
        "1 a1 1 ~t2_1 >= 1 ;\n"
        "1 ~a1 1 y0_1 >= 1 ;\n"
        "1 ~a1 1 y1_0 >= 1 ;\n"
        "1 ~a2 1 t4_1 >= 1 ;\n"
        "1 a2 1 ~t4_1 >= 1 ;\n"
        "1 ~a2 1 t2_0 1 t2_1 >= 1 ;\n"
        "1 ~t2_0 1 a2 >= 1 ;\n"
        "1 ~t2_1 1 a2 >= 1 ;\n"
        "1 ~t2_0 1 ~t2_1 >= 1 ;\n"
        "1 ~a3 >= 1 ;\n"
        "1 ~a4 1 t4_0 1 t4_1 >= 1 ;\n"
        "1 ~t4_0 1 a4 >= 1 ;\n"
        "1 ~t4_1 1 a4 >= 1 ;\n"
        "1 ~t4_0 1 ~t4_1 >= 1 ;\n");
}

TEST_CASE("the tagged formula has as many solutions as the statement's count")
{
    auto text = GENERATE(
        one_of_two,
        // overlapping children count twice
        statement("2", "(0 1 cube (0 0)) (1 1 cube (0 0)) (2 1 cube (0 1)) (3 3 or 0 1 2)", "3"),
        // a product of sums
        statement("3 2",
            "(0 1 cube (0 0)) (1 1 cube (0 2)) (2 2 or 0 1) (3 1 cube (1 0)) (4 1 cube (1 1)) (5 2 or 3 4) (6 4 and 2 5)", "6"),
        // a node shared by two ors, one of which overlaps it, and a dead node
        statement("2 2",
            "(0 1 cube (0 0) (1 1)) (1 1 cube (0 1) (1 0)) (2 2 or 0 1) (3 1 and) (4 3 or 0 2)", "4"),
        // an or below an and, beside a cube that overlaps it
        statement("2 2",
            "(0 1 cube (0 0)) (1 1 cube (0 1)) (2 2 or 0 1) (3 1 cube (1 1)) (4 2 and 2 3) "
            "(5 1 cube (0 0) (1 0)) (6 3 or 4 5)",
            "6"),
        // empty cubes and ands under an and
        statement("2", "(0 1 cube) (1 1 and) (2 1 cube (0 1)) (3 1 and 0 1 2)", "3"),
        // distinct values, all of them, a set of one, and one variable
        statement("4 4 4", "(0 24 distinct (0 1 2) (0 1 2 3))", "0"),
        statement("3 3 3", "(0 6 distinct (0 1 2) (0 1 2))", "0"),
        statement("4 4", "(0 3 distinct (0) (1 2 3)) (1 1 cube (1 0)) (2 3 and 0 1)", "2"),
        // under an or, beside a cube that overlaps it, and shared
        statement("3 3",
            "(0 2 distinct (0 1) (1 2)) (1 1 cube (0 1) (1 2)) (2 3 or 0 1) (3 6 distinct (0 1) (0 1 2)) "
            "(4 8 or 0 3) (5 11 or 2 4)",
            "5"));

    INFO(text);
    auto s = read(text);
    stringstream opb;
    write_tagged_opb(opb, s);
    stringstream expected;
    expected << s.nodes[s.root].count;
    CHECK(std::to_string(count_solutions(opb.str())) == expected.str());
}

TEST_CASE("counting statement writer writes what the reader reads")
{
    stringstream out;
    CountingStatementWriter writer{out, {3, 3}};
    auto a = writer.add_cube({{0, 0}});
    auto b = writer.add_cube({{0, 2}});
    auto c = writer.add_or({a, b});
    auto d = writer.add_cube({{1, 1}});
    auto e = writer.add_and({c, d});
    auto f = writer.add_cube({{0, 0}, {1, 1}});
    auto g = writer.add_or({e, f});
    writer.finish(g);
    CHECK(writer.number_of_nodes() == 7);
    CHECK(writer.count(g) == 3);

    auto s = read(out.str());
    CHECK(s.root == g);
    CHECK(s.nodes[g].count == 3);

    CHECK_THROWS_AS(writer.add_or({e, c}), CountingStatementError);
    CHECK_THROWS_AS(writer.add_and({a, a}), CountingStatementError);
}

namespace
{
    // Six maps of an edge 0-1 and a vertex 2 to four values: the edge to (0 1) or (1 0) and
    // the vertex to 2 or 3, or the edge to (2 3) and the vertex to 0 or 1.
    const string edge_and_vertex = statement("4 4 4",
        "(0 1 cube (0 0) (1 1)) (1 1 cube (0 1) (1 0)) (2 2 or 0 1) "
        "(3 1 cube (2 2)) (4 1 cube (2 3)) (5 2 or 3 4) (6 4 and 2 5) "
        "(7 1 cube (0 2) (1 3)) (8 1 cube (2 0)) (9 1 cube (2 1)) (10 2 or 8 9) (11 2 and 7 10) "
        "(12 6 or 6 11)",
        "12");
}

TEST_CASE("unranking a statement lists each of its objects once")
{
    auto s = read(edge_and_vertex);
    set<vector<int>> seen;
    for (long i = 0; i < 6; ++i) {
        auto values = unrank_counting_statement(s, loooong{i});
        REQUIRE(values.size() == 3);
        seen.insert(values);
    }
    CHECK(seen == set<vector<int>>{{0, 1, 2}, {0, 1, 3}, {1, 0, 2}, {1, 0, 3}, {2, 3, 0}, {2, 3, 1}});
    CHECK_THROWS_AS(unrank_counting_statement(s, loooong{6}), CountingStatementError);
    CHECK_THROWS_AS(unrank_counting_statement(s, loooong{-1}), CountingStatementError);
}

TEST_CASE("sampling a statement gives only its objects, and all of them")
{
    auto s = read(edge_and_vertex);
    CountingStatementSampler sample{s, 1};
    map<vector<int>, int> seen;
    for (int i = 0; i < 6000; ++i)
        ++seen[sample()];
    REQUIRE(seen.size() == 6);
    for (auto & [values, n] : seen) {
        CHECK(n > 800);
        CHECK(n < 1200);
    }
}

TEST_CASE("rendering a statement names products and factors tables")
{
    auto s = read(edge_and_vertex);
    stringstream out;
    render_counting_statement(out, s);
    CHECK(out.str() == "root [6] = one of\n"
                       "  | n2 * n5  [4]\n"
                       "  | (0->2 1->3) * 2 in {0 1}  [2]\n"
                       "n2 [2] = {0 1} = {0 1}  [2]\n"
                       "n5 [2] = 2 in {2 3}  [2]\n");
}

TEST_CASE("a distinct node's objects are every way of giving its variables distinct values")
{
    // vertex 2 at 3, and vertices 0 and 1 at distinct values of {0, 1, 2}; or vertex 2 at
    // 0, and vertices 0 and 1 at 1 and 2 either way round
    auto s = read(statement("4 4 4",
        "(0 6 distinct (0 1) (0 1 2)) (1 1 cube (2 3)) (2 6 and 0 1) "
        "(3 2 distinct (0 1) (1 2)) (4 1 cube (2 0)) (5 2 and 3 4) (6 8 or 2 5)",
        "6"));
    set<vector<int>> seen;
    for (long i = 0; i < 8; ++i)
        seen.insert(unrank_counting_statement(s, loooong{i}));
    CHECK(seen == set<vector<int>>{{0, 1, 3}, {0, 2, 3}, {1, 0, 3}, {1, 2, 3}, {2, 0, 3}, {2, 1, 3}, {1, 2, 0}, {2, 1, 0}});

    stringstream out;
    render_counting_statement(out, s);
    CHECK(out.str() == "root [8] = one of\n"
                       "  | (2->3) * {0 1} distinct in {0 1 2}  [6]\n"
                       "  | (2->0) * {0 1} = {1 2}  [2]\n");
}
