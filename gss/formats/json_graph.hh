#ifndef GLASGOW_SUBGRAPH_SOLVER_SOLVER_FORMATS_JSON_GRAPH_HH
#define GLASGOW_SUBGRAPH_SOLVER_SOLVER_FORMATS_JSON_GRAPH_HH 1

#include <gss/formats/graph_file_error.hh>
#include <gss/formats/input_graph.hh>

#include <iosfwd>
#include <string>

/**
 * The "gss-graph" JSON format. Unlike every other format we read, this one
 * declares its graph-level properties rather than leaving them to be inferred
 * from which syntax turned up in the file, and validates strictly: anything
 * ambiguous is an error rather than a guess.
 *
 * \code
 * {"format": "gss-graph", "version": 1, "directed": false,
 *  "vertices": 3, "edges": [[0, 1], [1, 2]]}
 * \endcode
 *
 * Vertices are either a count (that many anonymous vertices, addressed by index)
 * or an array whose entries are names, or objects with optional "name" and
 * "label". Edges are either the positional [from, to] / [from, to, label] — frozen
 * at two or three elements, so the compact form can never take on a new meaning —
 * or objects with "from", "to" and optional "label". Endpoints address a vertex by
 * index when they are integers and by name when they are strings, and a file has
 * to pick one. Keys beginning "x-" are ignored wherever they appear; any other
 * unrecognised key is an error.
 *
 * \throw GraphFileError for anything the format disallows, naming the offending
 *     key or index.
 */
auto read_json_graph(std::istream && infile, const std::string & filename) -> InputGraph;

/**
 * Write a graph in the "gss-graph" JSON format.
 *
 * The output is canonical: reading it back and writing it again gives the same
 * bytes, and every property the format carries (directedness, names, which label
 * kinds are present, loops) survives the trip. That round trip is what makes the
 * claim that the format is unambiguous testable, which is why the writer lives
 * alongside the reader rather than arriving later.
 */
auto write_json_graph(std::ostream & outfile, const InputGraph & graph) -> void;

#endif
