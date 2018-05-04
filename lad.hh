/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GUARD_LAD_HH
#define GUARD_LAD_HH 1

#include "graph.hh"
#include "graph_file_error.hh"

#include <string>

/**
 * Read a LAD format file into a Graph.
 *
 * \throw GraphFileError
 */
auto read_lad(const std::string & filename) -> Graph;

#endif
