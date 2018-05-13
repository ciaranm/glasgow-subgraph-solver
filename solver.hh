/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_SOLVER_HH
#define GLASGOW_SUBGRAPH_SOLVER_SOLVER_HH 1

#include "params.hh"
#include "result.hh"
#include "formats/input_graph.hh"

auto sequential_subgraph_isomorphism(const std::pair<InputGraph, InputGraph> & graphs, const Params & params) -> Result;

#endif
