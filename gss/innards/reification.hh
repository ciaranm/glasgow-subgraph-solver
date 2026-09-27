#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_GSS_INNARDS_REIFICATION_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_GSS_INNARDS_REIFICATION_HH 1

#include <gss/formats/input_graph.hh>

#include <vector>

namespace gss::innards
{
    /**
     * Where an edge-vertex came from: the endpoints of the edge it stands for, as
     * vertices of the reified graph (which are also the original vertex numbers). In a
     * graph reified as undirected, from and to are in no particular order.
     */
    struct EdgeVertexEndpoints
    {
        int from, to;
    };

    /**
     * A pattern and target with every edge turned into a vertex of its own, so that
     * parallel edges and edge costs become things the homomorphism solver already
     * handles: distinct vertices, and costs on target vertices.
     *
     * Each reified graph is a simple graph with vertex labels and no edge labels. Its
     * first vertices are the original vertices, with the same numbers, names and costs,
     * so a mapping of the reified graphs restricted to them is a mapping of the
     * originals. After them come the edge-vertices, one for each edge, joined to the
     * edge's endpoints: from -> e -> to if either graph is directed (so an undirected
     * edge of an undirected graph becomes two edge-vertices, one per direction, as it is
     * two arcs), and from - e - to otherwise.
     *
     * Labels are rewritten so that the two kinds of vertex can never be confused,
     * whatever text a user chose: an original vertex is labelled "v" followed by its
     * label, an edge-vertex "e" followed by its edge's label, and an edge-vertex for a
     * loop "l" followed by its label. So a pattern loop can only land on a target loop,
     * and a pattern edge only on a target edge. As in the model, a pattern with no vertex
     * labels (or no edge labels) ignores the target's, so they are left out of the
     * rewritten labels on both sides.
     *
     * An edge-vertex costs what its edge did, and an original vertex what it did, or 0
     * for a graph without costs of that kind.
     *
     * Only an injective, non-induced mapping of the reified graphs means the same thing
     * as a mapping of the originals: induced would talk about non-edges between
     * original vertices and edge-vertices, and a non-injective mapping may collapse a
     * pattern edge onto a target loop, which the loop labels forbid. solve_homomorphism_
     * problem refuses those combinations.
     */
    struct Reification
    {
        InputGraph pattern, target;
        int pattern_original_size, target_original_size;

        // Indexed by edge-vertex number minus the original size.
        std::vector<EdgeVertexEndpoints> pattern_edge_vertices, target_edge_vertices;

        bool directed;

        // Over every vertex of the reified target.
        std::vector<long long> target_costs;
    };

    /**
     * Does solving these graphs need reifying them first? True if either is a
     * multigraph, or if we are minimising cost over a target with edge costs.
     */
    auto needs_reification(const InputGraph & pattern, const InputGraph & target, bool minimising_cost) -> bool;

    auto reify(const InputGraph & pattern, const InputGraph & target) -> Reification;
}

#endif
