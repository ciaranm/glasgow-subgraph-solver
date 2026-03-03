#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_FORMATS_READ_FILE_FORMAT_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_FORMATS_READ_FILE_FORMAT_HH 1

#include <gss/formats/graph_file_error.hh>
#include <gss/formats/input_graph.hh>

#include <string>

/**
 * Detect a graph file format.
 *
 * \throw GraphFileError
 */
auto detect_file_format(std::istream & infile, const std::string & filename) -> std::string;

/**
 * Read in a file in the specified format ("auto" to try to auto-detect).
 *
 * \throw GraphFileError
 */
auto read_file_format(const std::string & format, const std::string & filename) -> InputGraph;

/**
 * Read in an already-open file in the specified format ("auto" requires seekable stream).
 *
 * \throw GraphFileError
 */
auto read_file_format(const std::string & format, const std::string & filename, std::istream & stream) -> InputGraph;

#endif
