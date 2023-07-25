#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_GRAPH_TRAITS_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_GRAPH_TRAITS_HH 1

#include <gss/formats/input_graph.hh>

namespace gss::innards
{
    auto is_simple_clique(const InputGraph & graph) -> bool;
}

#endif
