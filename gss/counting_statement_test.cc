#include <gss/counting_statement.hh>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <sstream>
#include <string>
#include <vector>

using namespace gss;

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
        return "((version 1) (domain " + domain + ") (nodes " + nodes + ") (root " + root + "))";
    }

    // one variable, three values, either of the first two
    const string one_of_two = statement("3", "(0 cube (0 0)) (1 cube (0 1)) (2 or 0 1)", "2");
}

TEST_CASE("counting statements read and count")
{
    auto s = read(one_of_two);
    CHECK(s.domain_sizes == vector<int>{3});
    CHECK(s.nodes.size() == 3);
    CHECK(count_counting_statement(s) == 2);
    CHECK(counting_side_conditions(s).empty());
}

TEST_CASE("counting statements multiply, add and subtract")
{
    // x0 in {0, 1}, x1 in {0, 1, 2}, minus the one where both are 0
    auto s = read(statement("3 3",
        "(0 cube (0 0)) (1 cube (0 1)) (2 or 0 1) "
        "(3 cube (1 0)) (4 cube (1 1)) (5 cube (1 2)) (6 or 3 4 5) "
        "(7 and 2 6) (8 cube (0 0) (1 0)) (9 minus 7 8)",
        "9"));
    auto counts = count_counting_nodes(s);
    CHECK(counts[2] == 2);
    CHECK(counts[6] == 3);
    CHECK(counts[7] == 6);
    CHECK(counts[9] == 5);
    CHECK(counting_side_conditions(s) == vector<int>{9});
}

TEST_CASE("counting statements allow empty cubes and ands")
{
    CHECK(count_counting_statement(read(statement("", "(0 cube)", "0"))) == 1);
    CHECK(count_counting_statement(read(statement("", "(0 and)", "0"))) == 1);
}

TEST_CASE("counting statements leave only non-syntactic side conditions")
{
    // an or of distinct cubes is discharged; an or with a repeated cube, or over
    // anything that is not a cube, is not
    auto s = read(statement("3",
        "(0 cube (0 0)) (1 cube (0 1)) (2 cube (0 1)) "
        "(3 or 0 1) (4 or 1 2) (5 or 0 0) (6 or 3)",
        "6"));
    CHECK(counting_side_conditions(s) == vector<int>{4, 5, 6});
}

TEST_CASE("counting statements reject a negative minus")
{
    auto s = read(statement("3", "(0 cube (0 0)) (1 cube (0 1)) (2 or 0 1) (3 minus 0 2)", "3"));
    CHECK_THROWS_AS(count_counting_statement(s), CountingStatementError);
}

TEST_CASE("counting statements reject anything malformed")
{
    auto bad = GENERATE(
        // lexical and grammatical
        string{""},
        string{"("},
        statement("3", "(0 cube (0 0)) (1 cube (0 1)) (2 or 0 1)", "2") + " extra",
        string{"((version 2) (domain 1 1) (nodes (0 cube (0 0))) (root 0))"},
        string{"((domain 1 1) (version 1) (nodes (0 cube (0 0))) (root 0))"},
        string{"((version 1) (domain 1 1) (nodes (0 cube (0 0))))"},
        statement("03", "(0 cube (0 0))", "0"),
        statement("-3", "(0 cube (0 0))", "0"),
        statement("+3", "(0 cube (0 0))", "0"),
        statement("3", "(0 cube (0 0)) (1 xor 0)", "1"),
        statement("3", "(0 cube (0 0) (0))", "0"),
        statement("3", "(0 cube (0 0 1))", "0"),
        statement("3", "(0 cube 0 0)", "0"),
        statement("3", "(0 or (0 0))", "0"),
        statement("99999999999", "(0 cube (0 0))", "0"),
        // ids and references
        statement("3", "(1 cube (0 0))", "1"),
        statement("3", "(0 cube (0 0)) (0 cube (0 1))", "0"),
        statement("3", "(0 or 0)", "0"),
        statement("3", "(0 or 1) (1 cube (0 0))", "0"),
        statement("3", "(0 cube (0 0))", "1"),
        // ranges
        statement("3", "(0 cube (1 0))", "0"),
        statement("3", "(0 cube (0 3))", "0"),
        // cubes
        statement("3 3", "(0 cube (1 0) (0 0))", "0"),
        statement("3 3", "(0 cube (0 0) (0 1))", "0"),
        // scopes
        statement("3 3", "(0 cube (0 0)) (1 cube (0 1) (1 1)) (2 and 0 1)", "2"),
        statement("3 3", "(0 cube (0 0) (1 0)) (1 cube (0 1)) (2 or 0 1)", "2"),
        statement("3 3", "(0 cube (0 0) (1 0)) (1 or)", "0"),
        statement("3 3", "(0 cube (0 0) (1 0)) (1 cube (0 1)) (2 minus 0 1)", "2"),
        statement("3 3", "(0 cube (0 0) (1 0)) (1 minus 0)", "1"),
        statement("3 3", "(0 cube (0 0) (1 0)) (1 minus 0 0 0)", "1"),
        // the root must cover every variable
        statement("3 3", "(0 cube (0 0))", "0"),
        statement("3 3", "(0 cube (0 0) (1 0)) (1 cube (0 0))", "1"),
        // each variable has its own domain
        statement("2 3", "(0 cube (0 2) (1 0))", "0"),
        statement("3 0", "(0 cube (0 0))", "0"));

    INFO(bad);
    CHECK_THROWS_AS(read(bad), CountingStatementError);
}

TEST_CASE("counting statements translate exactly as specified")
{
    auto s = read(statement("2 2",
        "(0 cube (0 0) (1 1)) (1 cube (0 1) (1 0)) (2 or 0 1) (3 cube (0 0) (1 1)) (4 minus 2 3) (5 and)",
        "4"));

    stringstream opb;
    write_counting_opb(opb, s);
    CHECK(opb.str() ==
        "preserved: y0_0 y0_1 y1_0 y1_1 ;\n"
        "1 y0_0 1 y0_1 >= 1 ;\n"
        "1 y0_0 1 y0_1 <= 1 ;\n"
        "1 y1_0 1 y1_1 >= 1 ;\n"
        "1 y1_0 1 y1_1 <= 1 ;\n"
        "2 ~g0 1 y0_0 1 y1_1 >= 2 ;\n"
        "1 g0 1 ~y0_0 1 ~y1_1 >= 1 ;\n"
        "2 ~g1 1 y0_1 1 y1_0 >= 2 ;\n"
        "1 g1 1 ~y0_1 1 ~y1_0 >= 1 ;\n"
        "1 ~g2 1 g0 1 g1 >= 1 ;\n"
        "2 g2 1 ~g0 1 ~g1 >= 2 ;\n"
        "2 ~g3 1 y0_0 1 y1_1 >= 2 ;\n"
        "1 g3 1 ~y0_0 1 ~y1_1 >= 1 ;\n"
        "1 ~g4 1 g2 >= 1 ;\n"
        "1 ~g4 1 ~g3 >= 1 ;\n"
        "1 g4 1 ~g2 1 g3 >= 1 ;\n"
        "1 g5 >= 1 ;\n"
        "1 g4 >= 1 ;\n");

    stringstream side;
    CHECK(write_side_condition_opb(side, s));
    auto side_text = side.str();
    auto definitions = opb.str().substr(opb.str().find('\n') + 1);
    definitions = definitions.substr(0, definitions.size() - string{"1 g4 >= 1 ;\n"}.size());
    CHECK(side_text ==
        definitions +
            "1 ~v4 1 g3 >= 1 ;\n"
            "1 ~v4 1 ~g2 >= 1 ;\n"
            "1 v4 >= 1 ;\n");
}

TEST_CASE("counting statements give each variable its own domain")
{
    // a variable with one value, and one with two, either of which is allowed
    auto s = read(statement("1 2", "(0 cube (0 0) (1 0)) (1 cube (0 0) (1 1)) (2 or 0 1)", "2"));
    CHECK(s.domain_sizes == vector<int>{1, 2});
    CHECK(count_counting_statement(s) == 2);

    stringstream opb;
    write_counting_opb(opb, s);
    CHECK(opb.str() ==
        "preserved: y0_0 y1_0 y1_1 ;\n"
        "1 y0_0 >= 1 ;\n"
        "1 y0_0 <= 1 ;\n"
        "1 y1_0 1 y1_1 >= 1 ;\n"
        "1 y1_0 1 y1_1 <= 1 ;\n"
        "2 ~g0 1 y0_0 1 y1_0 >= 2 ;\n"
        "1 g0 1 ~y0_0 1 ~y1_0 >= 1 ;\n"
        "2 ~g1 1 y0_0 1 y1_1 >= 2 ;\n"
        "1 g1 1 ~y0_0 1 ~y1_1 >= 1 ;\n"
        "1 ~g2 1 g0 1 g1 >= 1 ;\n"
        "2 g2 1 ~g0 1 ~g1 >= 2 ;\n"
        "1 g2 >= 1 ;\n");
}

TEST_CASE("counting statements write no side-condition instance when there is nothing to prove")
{
    stringstream side;
    CHECK(! write_side_condition_opb(side, read(one_of_two)));
    CHECK(side.str().empty());
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
    auto g = writer.add_minus(e, f);
    writer.finish(g);
    CHECK(writer.number_of_nodes() == 7);

    auto s = read(out.str());
    CHECK(s.root == g);
    CHECK(count_counting_statement(s) == 1);
}
