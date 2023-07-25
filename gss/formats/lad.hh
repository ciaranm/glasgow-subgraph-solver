#ifndef GLASGOW_SUBGRAPH_SOLVER_SOLVER_FORMATS_LAD_HH
#define GLASGOW_SUBGRAPH_SOLVER_SOLVER_FORMATS_LAD_HH 1

#include <gss/formats/graph_file_error.hh>
#include <gss/formats/input_graph.hh>

#include <iosfwd>
#include <string>

/**
 * Read a LAD format file into an InputGraph.
 *
 * \throw GraphFileError
 */
auto read_lad(std::ifstream && infile, const std::string & filename) -> InputGraph;

/**
 * Read a LAD format file into an InputGraph, treating edges as directed.
 *
 * \throw GraphFileError
 */
auto read_directed_lad(std::ifstream && infile, const std::string & filename) -> InputGraph;

/**
 * Read a Labelled LAD format file into an InputGraph.
 *
 * \throw GraphFileError
 */
auto read_labelled_lad(std::ifstream && infile, const std::string & filename) -> InputGraph;

/**
 * Read a Vertex-Labelled LAD format file into an InputGraph.
 *
 * \throw GraphFileError
 */
auto read_vertex_labelled_lad(std::ifstream && infile, const std::string & filename) -> InputGraph;

#endif
