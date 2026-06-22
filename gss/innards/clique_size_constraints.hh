#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_CLIQUE_SIZE_CONSTRAINTS_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_CLIQUE_SIZE_CONSTRAINTS_HH 1

#include <gss/homomorphism.hh>
#include <gss/innards/processed_graphs_data.hh>

#include <list>
#include <string>
#include <vector>

namespace gss::innards
{
    class HomomorphismProofs;

    // The clique-size filtering caches (the --cliques feature): for each graph pair, the
    // size of the largest clique around each pattern / target vertex (lazily computed),
    // the best-known bounds used to short-circuit those clique solves, and per-solve
    // timing / node-count statistics. Held by the model but operated on by the free
    // functions below; deliberately separate from the recoded graphs and the proof.
    struct CliqueSizeData
    {
        unsigned max_graphs_for_clique_size_constraints = 0;
        bool has_pattern_cliques_sizes = false;
        std::vector<std::vector<int>> pattern_cliques_sizes, target_cliques_sizes, pattern_cliques_best_knowns, target_cliques_best_knowns;
        std::vector<int> largest_pattern_clique;
        std::list<std::string> pattern_cliques_build_times, pattern_cliques_solve_times, pattern_cliques_solve_find_nodes, pattern_cliques_solve_prove_nodes;
        std::list<std::string> target_cliques_build_times, target_cliques_solve_times, target_cliques_solve_find_nodes, target_cliques_solve_prove_nodes;
    };

    // Size the caches for the configured number of graph pairs. No-op unless
    // params.clique_size_constraints is set.
    auto init_clique_size_data(CliqueSizeData & data, const HomomorphismParams & params,
        unsigned max_graphs, unsigned pattern_size, unsigned target_size) -> void;

    // Is mapping pattern vertex p to target vertex t permitted by the clique-size bound?
    // Lazily computes the clique sizes it needs; on a violation it emits the no-clique
    // proof through `proofs` (when non-null). Returns true (allowed) when clique-size
    // constraints are off.
    auto check_clique_compatibility(CliqueSizeData & data, const ProcessedGraphsData & graphs,
        unsigned max_graphs, unsigned pattern_size, unsigned target_size, const HomomorphismParams & params,
        HomomorphismProofs * proofs, int p, int t) -> bool;
}

#endif
