/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_SOLVER_FORMATS_CSV_HH
#define GLASGOW_SUBGRAPH_SOLVER_SOLVER_FORMATS_CSV_HH 1

#include "formats/input_graph.hh"
#include "formats/graph_file_error.hh"

#include <iosfwd>
#include <string>

/**
 * Read a CSV format file into an InputGraph.
 *
 * \throw GraphFileError
 */
auto read_csv(std::ifstream && infile, const std::string & filename) -> InputGraph;

auto read_csv_name(std::ifstream && infile, const std::string & filename, const std::string & name_map_filename) -> InputGraph;

#endif
