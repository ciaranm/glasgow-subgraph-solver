#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_MODEL_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_MODEL_HH 1

#include <gss/formats/input_graph.hh>
#include <gss/homomorphism.hh>
#include <gss/innards/homomorphism_domain.hh>
#include <gss/innards/proof.hh>
#include <gss/innards/svo_bitset.hh>

#include <array>
#include <cstdint>
#include <memory>
#include <vector>

namespace gss::innards
{
    class HomomorphismProofs;

    class HomomorphismModel
    {
    private:
        struct Imp;
        std::unique_ptr<Imp> _imp;

        auto _check_degree_compatibility(
            int p,
            int t,
            unsigned graphs_to_consider,
            std::vector<std::vector<std::vector<int>>> & patterns_ndss,
            std::vector<std::vector<std::optional<std::vector<int>>>> & targets_ndss,
            bool do_not_do_nds_yet) const -> bool;

        auto _check_loop_compatibility(int p, int t) const -> bool;

        auto _check_label_compatibility(int p, int t) const -> bool;

        auto _check_clique_compatibility(int p, int t) const -> bool;

        // Refresh the searcher's flat copy of active_graphs() (minus the original graph).
        auto _sync_active_supplemental_graphs() -> void;

    public:
        using PatternAdjacencyBitsType = uint8_t;

        const unsigned max_graphs;
        unsigned pattern_size, target_size;

        // The supplemental slots (g >= 1) the searcher still has to intersect with, and how
        // many there are -- active_graphs() without the original graph. A plain array member
        // rather than the vector behind it, because this is read once per
        // propagate_adjacency_constraints call: going through active_graphs() costs an
        // out-of-line call and two dependent loads (unique_ptr<Imp>, then the vector's heap
        // buffer) on a search that can run 65k nodes/s, which measurably outweighs the saving
        // on the targets where nothing is subsumed -- and those are the common case.
        // max_graphs is capped at 8 * sizeof(PatternAdjacencyBitsType) at construction.
        std::array<std::uint8_t, 8 * sizeof(PatternAdjacencyBitsType)> active_supplemental_graphs{};
        unsigned n_active_supplemental_graphs = 0;

        auto has_less_thans() const -> bool;
        auto has_occur_less_thans() const -> bool;
        std::vector<std::pair<unsigned, unsigned>> pattern_less_thans_in_convenient_order, target_occur_less_thans_in_convenient_order;

        HomomorphismModel(const InputGraph & target, const InputGraph & pattern, const HomomorphismParams & params,
            const std::shared_ptr<Proof> & proof, HomomorphismProofs * proofs);
        ~HomomorphismModel();

        auto pattern_vertex_for_proof(int v) const -> NamedVertex;
        auto target_vertex_for_proof(int v) const -> NamedVertex;

        // The solver-proofs middle layer (nullptr when not proving), so the searcher can
        // drive lazy supplemental materialisation on assignment / forward-check removal.
        [[nodiscard]] auto proofs() const -> HomomorphismProofs *;

        auto prepare() -> bool;

        // Build the supplemental graphs (exact-path, distance-3, …) into their slots and,
        // when proving, derive them. prepare() calls this itself unless staging is on, in
        // which case the staged solve driver calls it after a first bounded search round.
        // Precondition: prepare() has run and the original-graph self-loops are present.
        auto build_supplemental_graphs() -> void;

        // Tighten already-initialised domains using the supplemental graphs and NDS (the
        // filtering deferred past Stage 1 under staging), emitting the new prunings to the
        // proof. Returns false if a domain wipes out. Precondition: build_supplemental_graphs().
        auto tighten_domains_with_supplementals(std::vector<HomomorphismDomain> & domains) const -> bool;

        // The graph slots worth filtering with, ascending, always starting with the original
        // graph 0. A supplemental slot is dropped when an earlier slot in the same exact-path
        // run already subsumes it (identical target graph, and the exact-path graphs nest on
        // the pattern side) -- see build_supplemental_graphs. max_graphs, the bitset stride,
        // is unaffected: the slot still exists and is still built and proved, it is just never
        // re-tested. Meaningful only after build_supplemental_graphs; before it (Stage 1 under
        // staging) every slot is listed, which is harmless since only graph 0 is consulted then.
        [[nodiscard]] auto active_graphs() const -> const std::vector<unsigned> &;

        auto pattern_adjacency_bits(int p, int q) const -> PatternAdjacencyBitsType;
        auto pattern_graph_row(int g, int p) const -> const SVOBitset &;
        auto target_graph_row(int g, int t) const -> const SVOBitset &;

        auto forward_target_graph_row(int t) const -> const SVOBitset &;
        auto reverse_target_graph_row(int t) const -> const SVOBitset &;

        auto pattern_degree(int g, int p) const -> unsigned;
        auto target_degree(int g, int t) const -> unsigned;
        auto largest_target_degree() const -> unsigned;

        auto has_vertex_labels() const -> bool;
        auto has_edge_labels() const -> bool;
        auto directed() const -> bool;
        auto pattern_vertex_label(int p) const -> int;
        auto target_vertex_label(int p) const -> int;
        auto pattern_edge_label(int p, int q) const -> int;
        auto target_edge_label(int t, int u) const -> int;

        auto pattern_has_loop(int p) const -> bool;
        auto target_has_loop(int t) const -> bool;

        // When stage1 is true (staging), only the original graph is considered and NDS is
        // skipped -- the cheapest filtering, run before any supplemental graph is built.
        auto initialise_domains(std::vector<HomomorphismDomain> & domains, bool stage1 = false) const -> bool;

        auto add_extra_stats(std::list<std::string> &) const -> void;
    };
}

#endif
