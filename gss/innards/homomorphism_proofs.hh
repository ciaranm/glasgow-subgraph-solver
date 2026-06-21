#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_PROOFS_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_PROOFS_HH 1

#include <gss/formats/input_graph.hh>
#include <gss/homomorphism.hh>
#include <gss/innards/processed_graphs_data.hh>
#include <gss/innards/proof.hh>

#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <utility>
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
        std::shared_ptr<Proof> _proof;
        std::vector<std::string> _pattern_names, _target_names;

        // For each pattern edge (p,q) whose supplemental adjacency constraint was emitted
        // under subsumption elision, the slot whose (strongest) constraint we kept. Used to
        // re-derive an elided weaker constraint on demand (see ensure_supplemental_adjacency).
        std::map<std::pair<int, int>, unsigned> _kept_supplemental_slot;
        // Supplemental adjacency lines created transiently for a degree/NDS check, to delete
        // once it is done.
        std::vector<std::tuple<int, int, int, int>> _pending_transient_adjacencies;

        // Derive the exact-path (G^[gx2]) adjacency constraint for one head: "p maps to t
        // implies q maps to a vertex reachable from t by enough 2-walks". g is the slot the
        // constraint is recorded under; between_p_and_q are the (capped) common neighbours of
        // p and q; n_t / d_n_t are the neighbours of t in the original / exact-path graph; and
        // two_away_from_t pairs each vertex two hops from t with their common neighbours with
        // t. Emits through the generic Proof primitives over the shared adjacency cache.
        auto emit_exact_path_graph(int g, int p, int q, const std::vector<int> & between_p_and_q,
            int t, const std::vector<int> & n_t,
            const std::vector<std::pair<int, std::vector<int>>> & two_away_from_t,
            const std::vector<int> & d_n_t) -> void;

    public:
        HomomorphismProofs(const std::shared_ptr<Proof> & proof, const InputGraph & pattern, const InputGraph & target);

        [[nodiscard]] auto pattern_vertex(int v) const -> NamedVertex;
        [[nodiscard]] auto target_vertex(int v) const -> NamedVertex;

        // Emit the OPB model: a set of variables per CP variable, the (local-)injectivity
        // constraints, the adjacency constraints (and, for induced, the non-adjacency
        // constraints), the preserved set, then finalise. Operates on the raw input graphs
        // (it runs before the model is built), so it is the solver's first proof step.
        auto emit_model(const InputGraph & pattern, const InputGraph & target, const HomomorphismParams & params) -> void;

        // Derive the loop-cancelled forms of the adjacency constraints (a PBP derivation the
        // degree / supplemental-graph / distance-3 pols rely on). Deferred out of emit_model
        // and run by the search step, after the cheap concluding steps have had their chance:
        // an early conclusion (pattern-too-big, target-loop, clique) then never pays for these
        // derivations, while an instance that reaches search emits a byte-identical proof.
        auto derive_loop_fixed_adjacencies() -> void;

        // Prove the supplemental graphs that have just been built into the
        // ProcessedGraphsData. max_graphs is the bitset stride; number_of_exact_path_graphs
        // / slot identify which graph(s) within the data the derivation is about.
        // exact_path_index_and_slot pairs each surviving exact-path graph's index (the
        // path-count threshold) with the slot it occupies; exact_path_1_slot is where the
        // "at least one 2-path" graph lives (needed to justify all of them). Decoupling
        // the index from the slot lets inert-graph elimination renumber the supplementals.
        //
        // The supplemental adjacency constraints nest: for a fixed (p,t,q), the target set
        // shrinks as the path-count threshold grows (exact-path-k+1 subset of exact-path-k,
        // both subset of distance3), and a smaller-set constraint with the same head
        // syntactically subsumes a larger-set one. So we only emit the strongest (set-minimal)
        // constraint per head: prove_exact_path_graphs emits, for each pattern edge (p,q),
        // just the highest exact-path index that holds, and returns the (p,q) pairs it covered
        // so prove_distance3_graphs can skip them (distance3 is always the weakest). The
        // weaker constraints a per-graph degree/NDS check needs are created-then-deleted there
        // (see proof.cc); search-RUP can always fall back to the kept stronger constraint.
        // When elide_subsumed is false, every exact-path index's constraint is emitted (and
        // the returned covered set is empty, so distance3 also emits everything) -- the
        // pre-optimisation behaviour, kept behind --no-proof-supplemental-subsumption.
        auto prove_exact_path_graphs(const ProcessedGraphsData & graphs, unsigned max_graphs,
            const std::vector<std::pair<int, unsigned>> & exact_path_index_and_slot, unsigned exact_path_1_slot,
            bool elide_subsumed) -> std::set<std::pair<int, int>>;
        auto prove_distance3_graphs(const ProcessedGraphsData & graphs, unsigned max_graphs, unsigned slot,
            const std::set<std::pair<int, int>> & covered_by_exact_path) -> void;
        auto prove_extra_shape(const ProcessedGraphsData & graphs, unsigned max_graphs, unsigned slot) -> void;

        // A degree/NDS check needs the graph-g adjacency constraint for head (p->t, q). If
        // subsumption elision dropped it, re-derive it (a weakening of the kept stronger
        // constraint) so the check can cite it; record it so it can be deleted afterwards.
        // No-op if the constraint is still present (it was the kept one, or elision is off).
        auto ensure_supplemental_adjacency(const ProcessedGraphsData & graphs, unsigned max_graphs,
            int g, int p, int q, int t) -> void;
        // Delete every adjacency line ensure_supplemental_adjacency created since the last call.
        auto forget_transient_supplemental_adjacencies() -> void;

        // Prove that pattern vertex p cannot map to target tt, because p sits in a bigger
        // clique (in graph pair g) than tt does. Runs the clique solver with proof
        // extension over tt's neighbourhood to certify the bound.
        auto prove_no_clique(const ProcessedGraphsData & graphs, unsigned max_graphs, unsigned pattern_size,
            unsigned target_size, const HomomorphismParams & params, unsigned g, int p, int tt) -> void;
    };
}

#endif
