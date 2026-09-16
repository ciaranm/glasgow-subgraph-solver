#include <gss/innards/graph_traits.hh>

auto gss::innards::is_simple_clique(const InputGraph & graph) -> bool
{
    // "Simple" means: everything the clique solver's view of a graph leaves out. It has no
    // notion of a label, of a loop, or of an edge direction -- it builds one bitset row per
    // vertex from for_each_edge, which for a directed graph is asymmetric while every bound
    // and branching rule it has assumes symmetry. A complete digraph would otherwise pass
    // the edge count below (issue #93).
    if (graph.has_vertex_labels() || graph.has_edge_labels() || graph.loopy() || graph.directed())
        return false;

    return (graph.size() * (graph.size() - 1)) == graph.number_of_directed_edges();
}
