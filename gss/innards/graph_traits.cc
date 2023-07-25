#include <gss/innards/graph_traits.hh>

auto gss::innards::is_simple_clique(const InputGraph & graph) -> bool
{
    if (graph.has_vertex_labels() || graph.has_edge_labels() || graph.loopy())
        return false;

    return (graph.size() * (graph.size() - 1)) == graph.number_of_directed_edges();
}
