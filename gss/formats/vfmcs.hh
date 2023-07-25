#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_FORMATS_VFMCS_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_FORMATS_VFMCS_HH 1

#include <gss/formats/graph_file_error.hh>
#include <gss/formats/input_graph.hh>

#include <iosfwd>
#include <string>

auto read_unlabelled_undirected_vfmcs(std::ifstream && infile, const std::string & filename) -> InputGraph;

auto read_vertex_labelled_undirected_vfmcs(std::ifstream && infile, const std::string & filename) -> InputGraph;

auto read_vertex_labelled_directed_vfmcs(std::ifstream && infile, const std::string & filename) -> InputGraph;

#endif
