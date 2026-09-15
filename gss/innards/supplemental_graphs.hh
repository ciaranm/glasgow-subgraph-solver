#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_SUPPLEMENTAL_GRAPHS_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_SUPPLEMENTAL_GRAPHS_HH 1

#include <gss/formats/input_graph.hh>
#include <gss/homomorphism.hh>
#include <gss/innards/processed_graphs_data.hh>

namespace gss::innards
{
    // Build the supplemental graphs into a ProcessedGraphsData. These are pure
    // solver-core operations with no proof dependency (the proof derivations live in
    // HomomorphismProofs). Each is called once for the pattern side (pattern = true,
    // which also records the graph names) and once for the target side. idx is the next
    // free slot and is advanced past the graph(s) just built; max_graphs is the bitset
    // stride; size is the number of vertices on the side being built.
    auto build_exact_path_graphs(ProcessedGraphsData & graphs, unsigned size, unsigned & idx,
        unsigned max_graphs, unsigned number_of_exact_path_graphs, bool directed, bool at_most, bool pattern) -> void;

    auto build_distance3_graphs(ProcessedGraphsData & graphs, unsigned size, unsigned & idx,
        unsigned max_graphs, bool pattern) -> void;

    auto build_k4_graphs(ProcessedGraphsData & graphs, unsigned size, unsigned & idx,
        unsigned max_graphs, bool pattern) -> void;

    auto build_extra_shape(ProcessedGraphsData & graphs, unsigned size, unsigned & idx,
        unsigned max_graphs, const HomomorphismParams & params, InputGraph & shape, bool injective, int count, bool pattern) -> void;
}

#endif
