/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_SOLVER_FORMATS_BIGRAPH_HH
#define GLASGOW_SUBGRAPH_SOLVER_SOLVER_FORMATS_BIGRAPH_HH 1

#include "formats/input_graph.hh"
#include "formats/graph_file_error.hh"

#include <iosfwd>
#include <string>

/**
 * Read a pattern Bigraph format file into an InputGraph.
 *
 * \throw GraphFileError
 */
auto read_pattern_bigraph(std::istream && infile, const std::string & filename, bool directed) -> InputGraph;

/**
 * Read a target Bigraph format file into an InputGraph.
 *
 * \throw GraphFileError
 */
auto read_target_bigraph(std::istream && infile, const std::string & filename, bool directed) -> InputGraph;

#endif
