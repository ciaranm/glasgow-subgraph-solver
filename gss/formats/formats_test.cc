#include <gss/formats/csv.hh>
#include <gss/formats/dimacs.hh>
#include <gss/formats/graph_file_error.hh>
#include <gss/formats/input_graph.hh>
#include <gss/formats/lad.hh>
#include <gss/formats/vfmcs.hh>

#include <catch2/catch_test_macros.hpp>

#include <sstream>
#include <stdexcept>
#include <string>

using std::string;
using std::stringstream;

namespace
{
    // Look a vertex up by name; fails the test loudly rather than dereferencing a
    // disengaged optional if the name is missing.
    auto id(const InputGraph & g, const string & name) -> int
    {
        auto v = g.vertex_from_name(name);
        REQUIRE(v.has_value());
        return *v;
    }
}

// ---------------------------------------------------------------------------
// InputGraph
// ---------------------------------------------------------------------------

TEST_CASE("InputGraph: directedness is declared, not inferred")
{
    // A directed graph with no edges at all is still directed. Nothing about the
    // edges can tell you that, which is why the constructor has to be told.
    InputGraph empty_directed{3, false, false, true};
    CHECK(empty_directed.directed());
    CHECK(empty_directed.number_of_directed_edges() == 0);

    InputGraph undirected{3, false, false};
    CHECK_FALSE(undirected.directed());

    // Labelling an undirected edge does not change that.
    undirected.add_edge(0, 1, "red");
    CHECK_FALSE(undirected.directed());
    CHECK(undirected.adjacent(1, 0));
}

TEST_CASE("InputGraph: add_directed_edge needs the graph declared directed")
{
    InputGraph undirected{2, false, false};
    CHECK_THROWS_AS(undirected.add_directed_edge(0, 1, ""), std::logic_error);

    InputGraph directed{2, false, false, true};
    CHECK_NOTHROW(directed.add_directed_edge(0, 1, ""));
    CHECK(directed.adjacent(0, 1));
    CHECK_FALSE(directed.adjacent(1, 0));
}

TEST_CASE("InputGraph: vertex_has_name tells a set name from the index fallback")
{
    InputGraph g{2, false, false};
    CHECK_FALSE(g.vertex_has_name(0));
    CHECK(g.vertex_name(0) == "0"); // the fallback, not a name

    g.set_vertex_name(0, "alice");
    CHECK(g.vertex_has_name(0));
    CHECK(g.vertex_name(0) == "alice");

    // A vertex named for its own index is still a named vertex.
    g.set_vertex_name(1, "1");
    CHECK(g.vertex_has_name(1));
}

// ---------------------------------------------------------------------------
// CSV
// ---------------------------------------------------------------------------

TEST_CASE("read_csv: undirected edges")
{
    auto g = read_csv(stringstream{"a,b\nb,c\n"}, "g");
    CHECK(g.size() == 3);
    CHECK_FALSE(g.directed());
    CHECK_FALSE(g.has_vertex_labels());
    CHECK_FALSE(g.has_edge_labels());
    CHECK(g.adjacent(id(g, "a"), id(g, "b")));
    CHECK(g.adjacent(id(g, "b"), id(g, "a"))); // undirected: both directions
    CHECK_FALSE(g.adjacent(id(g, "a"), id(g, "c")));
    CHECK(g.number_of_directed_edges() == 4); // two undirected edges
}

TEST_CASE("read_csv: directed edges")
{
    auto g = read_csv(stringstream{"a>b\n"}, "g");
    CHECK(g.directed());
    CHECK(g.adjacent(id(g, "a"), id(g, "b")));
    CHECK_FALSE(g.adjacent(id(g, "b"), id(g, "a")));
    CHECK(g.number_of_directed_edges() == 1);
}

TEST_CASE("read_csv: vertex labels")
{
    auto g = read_csv(stringstream{"a,,red\nb,,blue\na,b\n"}, "g");
    CHECK(g.size() == 2);
    CHECK(g.has_vertex_labels());
    CHECK(g.vertex_label(id(g, "a")) == "red");
    CHECK(g.vertex_label(id(g, "b")) == "blue");
}

TEST_CASE("read_csv: partial vertex labelling is an error")
{
    // 'b' declares no label. Silently treating that as the empty-string label would
    // constrain it to unlabelled target vertices, so it has to be rejected instead.
    CHECK_THROWS_AS(read_csv(stringstream{"a,,X\nb,\na,b\n"}, "g"), GraphFileError);

    // Same thing, but 'b' is only ever mentioned by an edge line.
    CHECK_THROWS_AS(read_csv(stringstream{"a,,X\na,b\n"}, "g"), GraphFileError);
}

TEST_CASE("read_csv: edge labels")
{
    auto g = read_csv(stringstream{"a,b,red\n"}, "g");
    CHECK(g.has_edge_labels());
    CHECK(g.edge_label(id(g, "a"), id(g, "b")) == "red");
    CHECK(g.edge_label(id(g, "b"), id(g, "a")) == "red");
    // A labelled undirected edge is still undirected: adding a label must not flip
    // directed().
    CHECK_FALSE(g.directed());
    CHECK(g.number_of_directed_edges() == 2);
}

TEST_CASE("read_csv: partial edge labelling is an error")
{
    // The b--c edge declares no label, which is not the same thing as labelling it
    // with the empty string, so it can't be quietly read as that.
    CHECK_THROWS_AS(read_csv(stringstream{"a,b,red\nb,c\n"}, "g"), GraphFileError);

    // Also when the unlabelled edge comes first, and when the edges are directed.
    CHECK_THROWS_AS(read_csv(stringstream{"a,b\nb,c,red\n"}, "g"), GraphFileError);
    CHECK_THROWS_AS(read_csv(stringstream{"a>b,red\nb>c\n"}, "g"), GraphFileError);
}

TEST_CASE("read_csv: labelling one element type does not require labelling the other")
{
    // Vertex labels with unlabelled edges...
    auto vertices_only = read_csv(stringstream{"a,,X\nb,,Y\na,b\n"}, "g");
    CHECK(vertices_only.has_vertex_labels());
    CHECK_FALSE(vertices_only.has_edge_labels());
    CHECK(vertices_only.adjacent(id(vertices_only, "a"), id(vertices_only, "b")));

    // ...and edge labels with unlabelled vertices.
    auto edges_only = read_csv(stringstream{"a,b,red\n"}, "g");
    CHECK_FALSE(edges_only.has_vertex_labels());
    CHECK(edges_only.has_edge_labels());
}

TEST_CASE("read_csv: labelling an undirected edge changes nothing but the label")
{
    auto unlabelled = read_csv(stringstream{"a,b\nb,c\n"}, "g");
    auto labelled = read_csv(stringstream{"a,b,red\nb,c,red\n"}, "g");

    CHECK(labelled.size() == unlabelled.size());
    CHECK(labelled.directed() == unlabelled.directed());
    CHECK(labelled.loopy() == unlabelled.loopy());
    CHECK(labelled.number_of_directed_edges() == unlabelled.number_of_directed_edges());
    CHECK(labelled.has_edge_labels());
    CHECK_FALSE(unlabelled.has_edge_labels());
}

TEST_CASE("read_csv: directed edges with labels are still directed")
{
    auto g = read_csv(stringstream{"a>b,red\n"}, "g");
    CHECK(g.directed());
    CHECK(g.has_edge_labels());
    CHECK(g.edge_label(id(g, "a"), id(g, "b")) == "red");
    CHECK_FALSE(g.adjacent(id(g, "b"), id(g, "a")));
    CHECK(g.number_of_directed_edges() == 1);
}

TEST_CASE("read_csv: self loops")
{
    auto g = read_csv(stringstream{"a,a\na,b\n"}, "g");
    CHECK(g.loopy());
    CHECK(g.adjacent(id(g, "a"), id(g, "a")));
}

TEST_CASE("read_csv: a line without a delimiter is an error")
{
    CHECK_THROWS_AS(read_csv(stringstream{"a,b\nnodelimiter\n"}, "g"), GraphFileError);
}

TEST_CASE("read_csv: non-printable characters in a name are rejected")
{
    CHECK_THROWS_AS(read_csv(stringstream{string("a\x01z,b\n")}, "g"), GraphFileError);
}

// ---------------------------------------------------------------------------
// LAD
// ---------------------------------------------------------------------------

TEST_CASE("read_lad: undirected graph")
{
    // 3 vertices; vertex 0 adjacent to 1 and 2; 1 and 2 adjacent to 0.
    auto g = read_lad(stringstream{"3  2 1 2  1 0  1 0"}, "g");
    CHECK(g.size() == 3);
    CHECK_FALSE(g.directed());
    CHECK(g.adjacent(0, 1));
    CHECK(g.adjacent(1, 0));
    CHECK(g.adjacent(0, 2));
    CHECK_FALSE(g.adjacent(1, 2));
    CHECK(g.vertex_name(0) == "0");
    CHECK(g.number_of_directed_edges() == 4);
}

TEST_CASE("read_lad: trailing assignments rename vertices")
{
    auto g = read_lad(stringstream{"2  1 1  1 0  0=alice 1=bob"}, "g");
    CHECK(g.vertex_name(0) == "alice");
    CHECK(g.vertex_name(1) == "bob");
    CHECK(id(g, "alice") == 0);
}

TEST_CASE("read_directed_lad: edges are directed")
{
    auto g = read_directed_lad(stringstream{"2  1 1  0"}, "g");
    CHECK(g.directed());
    CHECK(g.adjacent(0, 1));
    CHECK_FALSE(g.adjacent(1, 0));
    CHECK(g.number_of_directed_edges() == 1);
}

TEST_CASE("read_vertex_labelled_lad: reads vertex labels")
{
    // size; then per vertex: <label> <degree> <neighbour>...
    auto g = read_vertex_labelled_lad(stringstream{"2  5 1 1  7 1 0"}, "g");
    CHECK(g.has_vertex_labels());
    CHECK_FALSE(g.has_edge_labels());
    CHECK(g.vertex_label(0) == "5");
    CHECK(g.vertex_label(1) == "7");
    CHECK(g.adjacent(0, 1));
}

TEST_CASE("read_labelled_lad: reads vertex and edge labels")
{
    // size; then per vertex: <vlabel> <degree> (<neighbour> <edgelabel>)...
    auto g = read_labelled_lad(stringstream{"2  5 1 1 9  7 0"}, "g");
    CHECK(g.has_vertex_labels());
    CHECK(g.has_edge_labels());
    CHECK(g.directed()); // this format is directed in its own right
    CHECK(g.vertex_label(0) == "5");
    CHECK(g.edge_label(0, 1) == "9");
}

TEST_CASE("read_vertex_labelled_lad: vertex labels do not make the graph directed")
{
    auto g = read_vertex_labelled_lad(stringstream{"2  5 1 1  7 1 0"}, "g");
    CHECK_FALSE(g.directed());
    CHECK(g.number_of_directed_edges() == 2);
}

TEST_CASE("read_lad: an out-of-bounds edge is an error")
{
    CHECK_THROWS_AS(read_lad(stringstream{"2  1 5  0"}, "g"), GraphFileError);
}

TEST_CASE("read_lad: trailing junk is an error")
{
    CHECK_THROWS_AS(read_lad(stringstream{"2  1 1  1 0  notanassignment"}, "g"), GraphFileError);
}

// ---------------------------------------------------------------------------
// DIMACS
// ---------------------------------------------------------------------------

TEST_CASE("read_dimacs: basic graph, 1-indexed")
{
    auto g = read_dimacs(stringstream{"p edge 3 2\ne 1 2\ne 2 3\n"}, "g");
    CHECK(g.size() == 3);
    CHECK(g.adjacent(0, 1)); // edge "1 2" -> 0-indexed 0-1
    CHECK(g.adjacent(1, 2));
    CHECK_FALSE(g.adjacent(0, 2));
    CHECK(g.vertex_name(0) == "1"); // names are 1-indexed
}

TEST_CASE("read_dimacs: comments are ignored")
{
    auto g = read_dimacs(stringstream{"c a comment\nc another\np edge 2 1\ne 1 2\n"}, "g");
    CHECK(g.size() == 2);
    CHECK(g.adjacent(0, 1));
}

TEST_CASE("read_dimacs: an edge before the problem line is out of bounds")
{
    CHECK_THROWS_AS(read_dimacs(stringstream{"e 1 2\n"}, "g"), GraphFileError);
}

TEST_CASE("read_dimacs: multiple problem lines are an error")
{
    CHECK_THROWS_AS(read_dimacs(stringstream{"p edge 2 0\np edge 3 0\n"}, "g"), GraphFileError);
}

TEST_CASE("read_dimacs: an unparseable line is an error")
{
    CHECK_THROWS_AS(read_dimacs(stringstream{"p edge 2 1\nx 1 2\n"}, "g"), GraphFileError);
}

// ---------------------------------------------------------------------------
// VFMCS (little-endian 16-bit binary)
// ---------------------------------------------------------------------------

namespace
{
    auto vfmcs_word(string & out, unsigned w) -> void
    {
        out.push_back(static_cast<char>(w & 0xff));
        out.push_back(static_cast<char>((w >> 8) & 0xff));
    }
}

TEST_CASE("read_unlabelled_undirected_vfmcs: basic graph")
{
    // 2 vertices, one undirected edge between them. Layout: size; one attribute
    // word per vertex; then per vertex an edge count followed by (target, label)
    // word pairs.
    string bytes;
    vfmcs_word(bytes, 2); // size
    vfmcs_word(bytes, 0); // vertex 0 attribute
    vfmcs_word(bytes, 0); // vertex 1 attribute
    vfmcs_word(bytes, 1); // vertex 0 edge count
    vfmcs_word(bytes, 1); // vertex 0 -> 1
    vfmcs_word(bytes, 0); // edge label (ignored)
    vfmcs_word(bytes, 0); // vertex 1 edge count

    auto g = read_unlabelled_undirected_vfmcs(stringstream{bytes}, "g");
    CHECK(g.size() == 2);
    CHECK_FALSE(g.directed());
    CHECK(g.adjacent(0, 1));
    CHECK(g.adjacent(1, 0));
    CHECK(g.number_of_directed_edges() == 2);
}

TEST_CASE("read_unlabelled_undirected_vfmcs: an out-of-bounds edge is an error")
{
    string bytes;
    vfmcs_word(bytes, 1); // size 1
    vfmcs_word(bytes, 0); // vertex 0 attribute
    vfmcs_word(bytes, 1); // vertex 0 edge count
    vfmcs_word(bytes, 5); // vertex 0 -> 5 (out of bounds)
    vfmcs_word(bytes, 0); // edge label

    CHECK_THROWS_AS(read_unlabelled_undirected_vfmcs(stringstream{bytes}, "g"), GraphFileError);
}
