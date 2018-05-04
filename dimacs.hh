/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GUARD_DIMACS_HH
#define GUARD_DIMACS_HH 1

#include "graph.hh"
#include "graph_file_error.hh"

#include <string>

/**
 * Read a DIMACS format file into a Graph.
 *
 * \throw GraphFileError
 */
auto read_dimacs(const std::string & filename) -> Graph;

#endif
