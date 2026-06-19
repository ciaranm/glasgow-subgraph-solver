#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_PROOFS_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_PROOFS_HH 1

#include <gss/formats/input_graph.hh>
#include <gss/innards/processed_graphs_data.hh>
#include <gss/innards/proof.hh>

#include <string>
#include <vector>

namespace gss::innards
{
    // The "solver-proofs" middle layer for the homomorphism solver. It owns the mapping
    // from internal vertex indices to their proof names, and the homomorphism-specific
    // proof derivations, reading the ProcessedGraphsData directly rather than having the
    // solver core marshal it across a narrow boundary. It is the solver's interface to
    // proof logging (solver core -> HomomorphismProofs -> generic Proof writer).
    //
    // For now its derivations still call the generic Proof methods, which take pre-built
    // NamedVertex vectors -- so that marshalling lives here now. Leaning the Proof API so
    // these can walk the bitsets and emit directly (deleting the marshalling) is a later
    // step; see dev_docs/preprocessor-refactor.md.
    class HomomorphismProofs
    {
    private:
        Proof & _proof;
        std::vector<std::string> _pattern_names, _target_names;

    public:
        HomomorphismProofs(Proof & proof, const InputGraph & pattern, const InputGraph & target);

        [[nodiscard]] auto pattern_vertex(int v) const -> NamedVertex;
        [[nodiscard]] auto target_vertex(int v) const -> NamedVertex;

        // Prove the supplemental graphs that have just been built into the
        // ProcessedGraphsData. max_graphs is the bitset stride; number_of_exact_path_graphs
        // / slot identify which graph(s) within the data the derivation is about.
        auto prove_exact_path_graphs(const ProcessedGraphsData & graphs, unsigned max_graphs, int number_of_exact_path_graphs) -> void;
        auto prove_distance3_graphs(const ProcessedGraphsData & graphs, unsigned max_graphs, unsigned slot) -> void;
        auto prove_extra_shape(const ProcessedGraphsData & graphs, unsigned max_graphs, unsigned slot) -> void;
    };
}

#endif
