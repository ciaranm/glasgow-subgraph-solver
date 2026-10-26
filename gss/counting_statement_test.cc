#include <gss/counting_statement.hh>

#include <catch2/catch_test_macros.hpp>

#include <sstream>
#include <string>

using namespace gss;

using std::string;
using std::stringstream;

namespace
{
    // Two variables over {p, q, r}: the four assignments where the first is p or the
    // second is q, as (p, anything) + (not p, q), less nothing.
    auto small_statement() -> CountingStatement
    {
        CountingStatement s;
        auto d = s.add_domain({"p", "q", "r"});
        auto u = s.add_variable("u", d);
        auto v = s.add_variable("v", d);
        auto up = s.add_cube({{u, 0}});
        auto any_v = s.add_sum({s.add_cube({{v, 0}}), s.add_cube({{v, 1}}), s.add_cube({{v, 2}})});
        auto first = s.add_product({up, any_v});
        auto second = s.add_sum({s.add_cube({{u, 1}, {v, 1}}), s.add_cube({{u, 2}, {v, 1}})});
        s.root = s.add_sum({first, second});
        return s;
    }

    auto round_trip(const CountingStatement & s) -> CountingStatement
    {
        stringstream text;
        write_counting_statement(text, s);
        return read_counting_statement(text);
    }

    auto count_text(const string & text) -> CountingStatementCount
    {
        stringstream in{text};
        return count_counting_statement(read_counting_statement(in));
    }
}

TEST_CASE("counting statements count")
{
    auto result = count_counting_statement(small_statement());
    CHECK(result.count == 5);
    REQUIRE(result.obligations.size() == 1);
    CHECK(result.obligations.front().kind == CountingObligation::Kind::Disjoint);
}

TEST_CASE("counting statements survive a round trip")
{
    auto s = small_statement();
    auto t = round_trip(s);
    CHECK(t.domains == s.domains);
    CHECK(t.variables == s.variables);
    REQUIRE(t.nodes.size() == s.nodes.size());
    for (unsigned n = 0; n < s.nodes.size(); ++n) {
        CHECK(t.nodes[n].kind == s.nodes[n].kind);
        CHECK(t.nodes[n].literals == s.nodes[n].literals);
        CHECK(t.nodes[n].variables == s.nodes[n].variables);
        CHECK(t.nodes[n].children == s.nodes[n].children);
    }
    CHECK(t.root == s.root);
    CHECK(count_counting_statement(t).count == 5);
}

TEST_CASE("differences and false count")
{
    auto result = count_text(R"(gss-counting-statement 1
* a comment
domain 0 2 p q
variable 0 u 0
cube 0 1 0 0
cube 1 1 0 1
sum 2 2 0 1
difference 3 2 0
false 4 1 0
sum 5 2 3 4
root 5
)");
    CHECK(result.count == 1);
    CHECK(result.obligations.size() == 2);
}

TEST_CASE("counting statements refuse what they can check is wrong")
{
    const string header = "gss-counting-statement 1\ndomain 0 2 p q\nvariable 0 u 0\nvariable 1 v 0\n";

    // each of these breaks one rule, which the message should mention
    auto check_refused = [&](const string & body, const string & why) {
        INFO(body);
        try {
            count_text(header + body);
            FAIL("accepted");
        }
        catch (const CountingStatementError & e) {
            CHECK(string{e.what()}.find(why) != string::npos);
        }
    };

    check_refused("cube 0 2 0 0 1 0\nproduct 1 2 0 0\nroot 1\n", "overlap");
    check_refused("cube 0 1 0 0\ncube 1 1 1 0\nsum 2 2 0 1\nroot 2\n", "different scopes");
    check_refused("cube 0 2 0 0 1 0\ncube 1 2 0 0 1 0\nsum 2 2 0 1\nroot 2\n", "same cube");
    check_refused("cube 0 2 0 0 1 0\ncube 1 2 0 1 1 0\nsum 2 2 0 1\ndifference 3 0 2\nroot 3\n", "negative");
    check_refused("cube 0 1 0 0\nroot 0\n", "every variable");
    check_refused("cube 0 2 0 0 1 0\nsum 1 1 2\nroot 1\n", "not defined before");
    check_refused("cube 0 2 0 0 0 1\nroot 0\n", "twice");
    check_refused("cube 0 2 0 2 1 0\nroot 0\n", "undefined value");
    check_refused("cube 0 2 0 0 2 0\nroot 0\n", "undefined variable");
    check_refused("sum 0 0\nroot 0\n", "no children");
    check_refused("cube 1 2 0 0 1 0\nroot 0\n", "consecutive");
    check_refused("cube 0 2 0 0 1 0\nroot 0\ncube 1 2 0 1 1 0\n", "follow the root");
    check_refused("cube 0 2 0 0 1 0\n", "missing root");
}
