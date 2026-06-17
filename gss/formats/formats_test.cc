#include <gss/formats/csv.hh>
#include <gss/formats/dimacs.hh>
#include <gss/formats/graph_file_error.hh>
#include <gss/formats/input_graph.hh>
#include <gss/formats/lad.hh>
#include <gss/formats/vfmcs.hh>

#include <catch2/catch_test_macros.hpp>

#include <sstream>
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

TEST_CASE("read_csv: edge labels")
{
    auto g = read_csv(stringstream{"a,b,red\n"}, "g");
    CHECK(g.has_edge_labels());
    CHECK(g.edge_label(id(g, "a"), id(g, "b")) == "red");
    CHECK(g.edge_label(id(g, "b"), id(g, "a")) == "red");
    // Quirk worth locking in: a labelled *undirected* CSV edge is built from two
    // add_directed_edge calls, so the graph reports itself as directed.
    CHECK(g.directed());
    CHECK(g.number_of_directed_edges() == 2);
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
    CHECK(g.directed());
    CHECK(g.vertex_label(0) == "5");
    CHECK(g.edge_label(0, 1) == "9");
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
