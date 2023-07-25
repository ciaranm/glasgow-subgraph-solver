#ifndef GLASGOW_SUBGRAPH_SOLVER_SOLVER_FORMATS_DIMACS_HH
#define GLASGOW_SUBGRAPH_SOLVER_SOLVER_FORMATS_DIMACS_HH 1

#include <gss/formats/graph_file_error.hh>
#include <gss/formats/input_graph.hh>

#include <iosfwd>
#include <string>

/**
 * Read a DIMACS format file into an InputGraph.
 *
 * \throw GraphFileError
 */
auto read_dimacs(std::ifstream && infile, const std::string & filename) -> InputGraph;

#endif
