#include <gss/formats/csv.hh>
#include <gss/formats/graph_file_error.hh>
#include <gss/formats/input_graph.hh>
#include <gss/formats/json_graph.hh>
#include <gss/formats/lad.hh>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <sstream>
#include <string>
#include <string_view>

using std::string;
using std::string_view;
using std::stringstream;

namespace
{
    auto parse(const string & text) -> InputGraph
    {
        return read_json_graph(stringstream{text}, "g");
    }

    auto write(const InputGraph & g) -> string
    {
        stringstream out;
        write_json_graph(out, g);
        return out.str();
    }

    auto id(const InputGraph & g, const string & name) -> int
    {
        auto v = g.vertex_from_name(name);
        REQUIRE(v.has_value());
        return *v;
    }

    // Everything the format is supposed to carry. Two graphs that agree here are
    // the same graph as far as anything downstream can tell.
    //
    // The label-presence flags are reported only for a graph that has elements of
    // that kind: "labelled, but with no vertices" and "labelled, but with no edges"
    // are states the flags can be left in by other readers, and neither is
    // expressible in a file or observable in a solve, there being nothing for the
    // labels to be attached to.
    auto describe(const InputGraph & g) -> string
    {
        stringstream s;
        s << "size=" << g.size()
          << " directed=" << g.directed()
          << " loopy=" << g.loopy()
          << " directed_edges=" << g.number_of_directed_edges();
        if (g.size() != 0)
            s << " vertex_labels=" << g.has_vertex_labels();
        if (g.number_of_directed_edges() != 0)
            s << " edge_labels=" << g.has_edge_labels();
        s << "\n";

        for (int v = 0; v < g.size(); ++v) {
            s << "  vertex " << v << " named=" << g.vertex_has_name(v) << " name=" << g.vertex_name(v);
            if (g.size() != 0 && g.has_vertex_labels())
                s << " label=" << g.vertex_label(v);
            s << "\n";
        }

        g.for_each_edge([&](int f, int t, string_view l) {
            s << "  edge " << f << " -> " << t << " label=" << l << "\n";
        });

        return s.str();
    }

    // The claim the format makes about itself: nothing is lost, and the writer's
    // output is canonical, so a second trip changes neither the graph nor the bytes.
    auto check_round_trip(const InputGraph & g) -> void
    {
        auto once = write(g);
        auto reparsed = parse(once);
        CHECK(describe(reparsed) == describe(g));
        CHECK(write(reparsed) == once);
    }
}

// ---------------------------------------------------------------------------
// Reading: the forms the format accepts
// ---------------------------------------------------------------------------

TEST_CASE("read_json_graph: minimal graph, vertices as a count")
{
    auto g = parse(R"({"format":"gss-graph","version":1,"directed":false,"vertices":3,"edges":[[0,1],[1,2]]})");
    CHECK(g.size() == 3);
    CHECK_FALSE(g.directed());
    CHECK_FALSE(g.loopy());
    CHECK_FALSE(g.has_vertex_labels());
    CHECK_FALSE(g.has_edge_labels());
    CHECK(g.adjacent(0, 1));
    CHECK(g.adjacent(1, 0)); // undirected: both ways
    CHECK_FALSE(g.adjacent(0, 2));
    CHECK(g.number_of_directed_edges() == 4);
    CHECK_FALSE(g.vertex_has_name(0));
}

TEST_CASE("read_json_graph: a directed graph keeps its direction")
{
    auto g = parse(R"({"format":"gss-graph","version":1,"directed":true,"vertices":2,"edges":[[0,1]]})");
    CHECK(g.directed());
    CHECK(g.adjacent(0, 1));
    CHECK_FALSE(g.adjacent(1, 0));
    CHECK(g.number_of_directed_edges() == 1);
}

TEST_CASE("read_json_graph: directed is declared, so it survives having no edges")
{
    // The property no other format of ours can express: there is nothing in an
    // empty edge list to infer directedness from.
    auto g = parse(R"({"format":"gss-graph","version":1,"directed":true,"vertices":3,"edges":[]})");
    CHECK(g.directed());
    CHECK(g.size() == 3);
    CHECK(g.number_of_directed_edges() == 0);
    check_round_trip(g);
}

TEST_CASE("read_json_graph: vertex names, as strings and as objects")
{
    auto named = parse(R"({"format":"gss-graph","version":1,"directed":false,
                           "vertices":["a","b"],"edges":[["a","b"]]})");
    CHECK(named.size() == 2);
    CHECK(named.vertex_has_name(0));
    CHECK(named.vertex_name(0) == "a");
    CHECK(named.adjacent(id(named, "a"), id(named, "b")));

    // ["a","b"] is exact sugar for [{"name":"a"},{"name":"b"}].
    auto verbose = parse(R"({"format":"gss-graph","version":1,"directed":false,
                             "vertices":[{"name":"a"},{"name":"b"}],"edges":[["a","b"]]})");
    CHECK(describe(verbose) == describe(named));
}

TEST_CASE("read_json_graph: vertex and edge labels")
{
    auto g = parse(R"({"format":"gss-graph","version":1,"directed":false,
                       "vertices":[{"name":"c1","label":"C"},{"name":"n1","label":"N"}],
                       "edges":[["c1","n1","single"]]})");
    CHECK(g.has_vertex_labels());
    CHECK(g.has_edge_labels());
    CHECK(g.vertex_label(id(g, "c1")) == "C");
    CHECK(g.vertex_label(id(g, "n1")) == "N");
    CHECK(g.edge_label(id(g, "c1"), id(g, "n1")) == "single");
    CHECK(g.edge_label(id(g, "n1"), id(g, "c1")) == "single");
    CHECK_FALSE(g.directed()); // a labelled undirected edge is still undirected
    check_round_trip(g);
}

TEST_CASE("read_json_graph: the object form of an edge")
{
    auto positional = parse(R"({"format":"gss-graph","version":1,"directed":false,
                                "vertices":2,"edges":[[0,1,"red"]]})");
    auto object = parse(R"({"format":"gss-graph","version":1,"directed":false,
                            "vertices":2,"edges":[{"from":0,"to":1,"label":"red"}]})");
    CHECK(describe(object) == describe(positional));

    // multiplicity 1 is what every ordinary edge has, so it is accepted.
    auto with_multiplicity = parse(R"({"format":"gss-graph","version":1,"directed":false,
                                       "vertices":2,"edges":[{"from":0,"to":1,"label":"red","multiplicity":1}]})");
    CHECK(describe(with_multiplicity) == describe(positional));
}

TEST_CASE("read_json_graph: a self-loop is exactly one edge")
{
    auto undirected = parse(R"({"format":"gss-graph","version":1,"directed":false,"vertices":2,"edges":[[0,0]]})");
    CHECK(undirected.loopy());
    CHECK(undirected.adjacent(0, 0));
    CHECK(undirected.number_of_directed_edges() == 1); // never implicitly doubled
    check_round_trip(undirected);

    auto directed = parse(R"({"format":"gss-graph","version":1,"directed":true,"vertices":2,"edges":[[0,0]]})");
    CHECK(directed.loopy());
    CHECK(directed.number_of_directed_edges() == 1);
    check_round_trip(directed);
}

TEST_CASE("read_json_graph: an isolated vertex needs no special syntax")
{
    auto g = parse(R"({"format":"gss-graph","version":1,"directed":false,
                       "vertices":["a","b","lonely"],"edges":[["a","b"]]})");
    CHECK(g.size() == 3);
    CHECK(g.degree(id(g, "lonely")) == 0);
    check_round_trip(g);
}

TEST_CASE("read_json_graph: the empty string is a real label")
{
    auto g = parse(R"({"format":"gss-graph","version":1,"directed":false,
                       "vertices":[{"label":""},{"label":""}],"edges":[[0,1,""]]})");
    CHECK(g.has_vertex_labels());
    CHECK(g.has_edge_labels());
    CHECK(g.vertex_label(0) == "");
    CHECK(g.edge_label(0, 1) == "");
    check_round_trip(g);
}

TEST_CASE("read_json_graph: named vertices may still be addressed by index")
{
    auto g = parse(R"({"format":"gss-graph","version":1,"directed":false,
                       "vertices":["a","b"],"edges":[[0,1]]})");
    CHECK(g.adjacent(id(g, "a"), id(g, "b")));
}

TEST_CASE("read_json_graph: a vertex may be named for a number without being an index")
{
    // The trap that costs every text format of ours: here the name "1" and the
    // index 1 are different tokens and cannot be confused.
    auto g = parse(R"({"format":"gss-graph","version":1,"directed":false,
                       "vertices":["1","0"],"edges":[["1","0"]]})");
    CHECK(id(g, "1") == 0);
    CHECK(id(g, "0") == 1);
    check_round_trip(g);
}

TEST_CASE("read_json_graph: x- keys are ignored wherever they appear")
{
    auto g = parse(R"({"format":"gss-graph","version":1,"directed":false,"x-provenance":"notes",
                       "vertices":[{"name":"a","x-colour":"blue"},{"name":"b"}],
                       "edges":[{"from":"a","to":"b","x-weight":2.5}]})");
    CHECK(g.size() == 2);
    CHECK(g.adjacent(id(g, "a"), id(g, "b")));
    CHECK_FALSE(g.has_vertex_labels());
}

// ---------------------------------------------------------------------------
// Reading: every rule that must fail, and be told why
// ---------------------------------------------------------------------------

namespace
{
    // Checks that the file is refused and that the message says what was wrong,
    // since a rejection nobody can act on is barely better than a wrong answer.
    auto check_rejected(const string & text, const string & expected_in_message) -> void
    {
        try {
            parse(text);
            FAIL("expected " + text + " to be rejected");
        }
        catch (const GraphFileError & e) {
            string message{e.what()};
            CHECK_THAT(message, Catch::Matchers::ContainsSubstring(expected_in_message));
        }
    }
}

TEST_CASE("read_json_graph: rule 1, format and version are checked")
{
    check_rejected(R"({"version":1,"directed":false,"vertices":0,"edges":[]})", "missing required key \"format\"");
    check_rejected(R"({"format":"something-else","version":1,"directed":false,"vertices":0,"edges":[]})", "only knows \"gss-graph\"");
    check_rejected(R"({"format":"gss-graph","directed":false,"vertices":0,"edges":[]})", "missing required key \"version\"");
    check_rejected(R"({"format":"gss-graph","version":"1","directed":false,"vertices":0,"edges":[]})", "\"version\" must be an integer");
    check_rejected(R"({"format":"gss-graph","version":0,"directed":false,"vertices":0,"edges":[]})", "not a version");
    // A newer file names both versions, so the reader that is too old says so.
    check_rejected(R"({"format":"gss-graph","version":2,"directed":false,"vertices":0,"edges":[]})", "version 2");
    check_rejected(R"({"format":"gss-graph","version":2,"directed":false,"vertices":0,"edges":[]})", "up to version 1");
}

TEST_CASE("read_json_graph: rule 2, directed is never inferred")
{
    check_rejected(R"({"format":"gss-graph","version":1,"vertices":2,"edges":[[0,1]]})", "missing required key \"directed\"");
    check_rejected(R"({"format":"gss-graph","version":1,"directed":"no","vertices":2,"edges":[[0,1]]})", "must be true or false");
}

TEST_CASE("read_json_graph: rule 3, a repeated edge is an error")
{
    check_rejected(R"({"format":"gss-graph","version":1,"directed":false,"vertices":2,"edges":[[0,1],[0,1]]})",
        "repeats the edge already given as edge 0");

    // In an undirected graph these are the same edge, not two.
    check_rejected(R"({"format":"gss-graph","version":1,"directed":false,"vertices":2,"edges":[[0,1],[1,0]]})",
        "[u, v] and [v, u] are the same edge");

    // In a directed graph they are genuinely different.
    auto both_ways = parse(R"({"format":"gss-graph","version":1,"directed":true,"vertices":2,"edges":[[0,1],[1,0]]})");
    CHECK(both_ways.number_of_directed_edges() == 2);

    check_rejected(R"({"format":"gss-graph","version":1,"directed":true,"vertices":2,"edges":[[0,1],[0,1]]})",
        "repeats the edge already given as edge 0");
}

TEST_CASE("read_json_graph: rule 5, labels are all or nothing per element type")
{
    check_rejected(R"({"format":"gss-graph","version":1,"directed":false,
                       "vertices":[{"name":"a","label":"X"},{"name":"b"}],"edges":[["a","b"]]})",
        "1 of 2 vertices carry a \"label\"");

    check_rejected(R"({"format":"gss-graph","version":1,"directed":false,
                       "vertices":3,"edges":[[0,1,"red"],[1,2]]})",
        "1 of 2 edges carry a \"label\"");

    // Labelling the vertices does not require labelling the edges, or the reverse.
    CHECK_NOTHROW(parse(R"({"format":"gss-graph","version":1,"directed":false,
                            "vertices":[{"label":"X"},{"label":"Y"}],"edges":[[0,1]]})"));
    CHECK_NOTHROW(parse(R"({"format":"gss-graph","version":1,"directed":false,
                            "vertices":2,"edges":[[0,1,"red"]]})"));
}

TEST_CASE("read_json_graph: rule 6, one addressing mode per file")
{
    check_rejected(R"({"format":"gss-graph","version":1,"directed":false,
                       "vertices":["a","b"],"edges":[["a","b"],[0,1]]})",
        "a file must address vertices one way throughout");

    check_rejected(R"({"format":"gss-graph","version":1,"directed":false,
                       "vertices":["a","b"],"edges":[["a","nope"]]})",
        "which is not in \"vertices\"");

    // Anonymous vertices have no names to address them by.
    check_rejected(R"({"format":"gss-graph","version":1,"directed":false,"vertices":2,"edges":[["a","b"]]})",
        "which is not in \"vertices\"");

    check_rejected(R"({"format":"gss-graph","version":1,"directed":false,"vertices":2,"edges":[[0,5]]})",
        "but the graph has 2 vertices");
    check_rejected(R"({"format":"gss-graph","version":1,"directed":false,"vertices":2,"edges":[[0,-1]]})",
        "but the graph has 2 vertices");
    check_rejected(R"({"format":"gss-graph","version":1,"directed":false,"vertices":2,"edges":[[0,1.5]]})",
        "must be an integer index or a string name");
}

TEST_CASE("read_json_graph: rule 7, an unrecognised key is an error")
{
    check_rejected(R"({"format":"gss-graph","version":1,"directed":false,"vertices":0,"edges":[],"wieght":1})",
        "unrecognised key \"wieght\"");
    check_rejected(R"({"format":"gss-graph","version":1,"directed":false,
                       "vertices":[{"name":"a","colour":"red"}],"edges":[]})",
        "unrecognised key \"colour\"");
    check_rejected(R"({"format":"gss-graph","version":1,"directed":false,
                       "vertices":2,"edges":[{"from":0,"to":1,"weight":3}]})",
        "unrecognised key \"weight\"");
}

TEST_CASE("read_json_graph: rule 8, vertex names are unique")
{
    check_rejected(R"({"format":"gss-graph","version":1,"directed":false,"vertices":["a","a"],"edges":[]})",
        "two vertices are both named \"a\"");
}

TEST_CASE("read_json_graph: the positional edge form is frozen at two or three elements")
{
    check_rejected(R"({"format":"gss-graph","version":1,"directed":false,"vertices":2,"edges":[[0,1,"red","extra"]]})",
        "has to use the object form");
    check_rejected(R"({"format":"gss-graph","version":1,"directed":false,"vertices":2,"edges":[[0]]})",
        "[from, to] or [from, to, label]");
}

TEST_CASE("read_json_graph: multi-edges are expressible but refused by this build")
{
    check_rejected(R"({"format":"gss-graph","version":1,"directed":false,"multigraph":true,"vertices":2,"edges":[[0,1]]})",
        "not supported by this build");
    check_rejected(R"({"format":"gss-graph","version":1,"directed":false,
                       "vertices":2,"edges":[{"from":0,"to":1,"multiplicity":3}]})",
        "not supported by this build");
    check_rejected(R"({"format":"gss-graph","version":1,"directed":false,
                       "vertices":2,"edges":[{"from":0,"to":1,"multiplicity":0}]})",
        "at least 1");
}

TEST_CASE("read_json_graph: structural rubbish is refused")
{
    check_rejected("[1,2,3]", "top level must be a JSON object");
    check_rejected("{not json", "could not be parsed as JSON");
    // Anything after the document is an error rather than being quietly ignored.
    check_rejected(R"({"format":"gss-graph","version":1,"directed":false,"vertices":0,"edges":[]} trailing)",
        "could not be parsed as JSON");
    check_rejected(R"({"format":"gss-graph","version":1,"directed":false,"edges":[]})",
        "missing required key \"vertices\"");
    check_rejected(R"({"format":"gss-graph","version":1,"directed":false,"vertices":0})",
        "missing required key \"edges\"");
    check_rejected(R"({"format":"gss-graph","version":1,"directed":false,"vertices":-1,"edges":[]})",
        "negative");
    check_rejected(R"({"format":"gss-graph","version":1,"directed":false,"vertices":{},"edges":[]})",
        "must be a count or an array");
    check_rejected(R"({"format":"gss-graph","version":1,"directed":false,"vertices":2,"edges":{}})",
        "\"edges\" must be an array");
    check_rejected(R"({"format":"gss-graph","version":1,"directed":false,"vertices":2,"edges":[7]})",
        "must be an array or an object");
    check_rejected(R"({"format":"gss-graph","version":1,"directed":false,"vertices":[7],"edges":[]})",
        "must be a name or an object");
    check_rejected(R"({"format":"gss-graph","version":1,"directed":false,"vertices":[{"name":7}],"edges":[]})",
        "must be a string");
    check_rejected(R"({"format":"gss-graph","version":1,"directed":false,"vertices":2,"edges":[[0,1,7]]})",
        "must be a string");
}

// ---------------------------------------------------------------------------
// Writing and round-tripping
// ---------------------------------------------------------------------------

TEST_CASE("write_json_graph: a graph with nothing to say about its vertices writes a count")
{
    auto g = parse(R"({"format":"gss-graph","version":1,"directed":false,"vertices":3,"edges":[[0,1],[1,2]]})");
    auto text = write(g);
    CHECK_THAT(text, Catch::Matchers::ContainsSubstring("\"vertices\": 3"));
    CHECK_THAT(text, Catch::Matchers::ContainsSubstring("\"directed\": false"));
    // An undirected edge is written once, not once per direction.
    CHECK_THAT(text, Catch::Matchers::ContainsSubstring("[0,1]"));
    CHECK_THAT(text, ! Catch::Matchers::ContainsSubstring("[1,0]"));
    check_round_trip(g);
}

TEST_CASE("write_json_graph: names are written, and used to address the edges")
{
    auto g = parse(R"({"format":"gss-graph","version":1,"directed":true,"vertices":["a","b"],"edges":[["a","b"]]})");
    auto text = write(g);
    CHECK_THAT(text, Catch::Matchers::ContainsSubstring(R"(["a","b"])"));
    CHECK_THAT(text, Catch::Matchers::ContainsSubstring("\"directed\": true"));
    check_round_trip(g);
}

TEST_CASE("write_json_graph: a partly named graph keeps which vertices were named")
{
    // Every other reader we have names either all of its vertices or none, so this
    // graph can only come from JSON -- and the writer has to represent it without
    // inventing names for the rest.
    auto g = parse(R"({"format":"gss-graph","version":1,"directed":false,
                       "vertices":[{"name":"alice"},{}],"edges":[[0,1]]})");
    REQUIRE(g.vertex_has_name(0));
    REQUIRE_FALSE(g.vertex_has_name(1));

    auto text = write(g);
    CHECK_THAT(text, Catch::Matchers::ContainsSubstring(R"({"name":"alice"})"));
    CHECK_THAT(text, Catch::Matchers::ContainsSubstring("{}")); // the unnamed one
    check_round_trip(g);
}

TEST_CASE("write_json_graph: round trip is a fixed point for graphs from other formats")
{
    SECTION("undirected CSV with edge labels")
    {
        check_round_trip(read_csv(stringstream{"a,b,red\nb,c,blue\n"}, "g"));
    }

    SECTION("directed CSV")
    {
        check_round_trip(read_csv(stringstream{"a>b\nb>c\n"}, "g"));
    }

    SECTION("CSV with vertex labels and a self-loop")
    {
        check_round_trip(read_csv(stringstream{"a,,X\nb,,Y\na,b\na,a\n"}, "g"));
    }

    SECTION("undirected LAD")
    {
        check_round_trip(read_lad(stringstream{"3  2 1 2  1 0  1 0"}, "g"));
    }

    SECTION("vertex labelled LAD")
    {
        check_round_trip(read_vertex_labelled_lad(stringstream{"2  5 1 1  7 1 0"}, "g"));
    }

    SECTION("labelled, directed LAD")
    {
        check_round_trip(read_labelled_lad(stringstream{"2  5 1 1 9  7 0"}, "g"));
    }
}
