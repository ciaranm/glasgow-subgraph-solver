#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_PROCESSED_GRAPHS_DATA_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_PROCESSED_GRAPHS_DATA_HH 1

#include <gss/innards/svo_bitset.hh>

#include <cstdint>
#include <list>
#include <string>
#include <vector>

namespace gss::innards
{
    // The two input graphs recoded as bitsets, plus everything derived from them that the
    // model's graph accessors expose: the per-graph adjacency rows (the original graph and
    // each supplemental graph, strided by the model's max_graphs), the degrees, loops,
    // labels, the compressed pattern adjacencies, and the names of the supplemental
    // graphs. This is the shared "processed graphs" that the supplemental builders, the
    // domain filters and the searcher all read.
    //
    // It deliberately holds no proof object, no vertex naming and no clique-size caches --
    // those stay private to the model (see dev_docs/preprocessor-refactor.md). The
    // dimensions (max_graphs / pattern_size / target_size) currently live on the model and
    // are passed in where the stride is needed; folding them in here rides the later
    // searcher-side restructuring.
    struct ProcessedGraphsData
    {
        std::vector<std::uint8_t> pattern_adjacencies_bits;
        std::vector<SVOBitset> pattern_graph_rows;
        std::vector<SVOBitset> target_graph_rows, forward_target_graph_rows, reverse_target_graph_rows;

        std::vector<std::vector<int>> patterns_degrees, targets_degrees;
        int largest_target_degree = 0;

        std::vector<int> pattern_vertex_labels, target_vertex_labels, pattern_edge_labels, target_edge_labels;
        std::vector<int> pattern_loops, target_loops;
        bool has_loops = false;
        bool directed = false;

        std::list<std::string> supplemental_graph_names;
    };
}

#endif
