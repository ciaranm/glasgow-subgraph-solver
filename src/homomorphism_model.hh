/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_MODEL_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_MODEL_HH 1

#include "formats/input_graph.hh"
#include "svo_bitset.hh"
#include "homomorphism.hh"
#include "homomorphism_domain.hh"
#include "proof.hh"

#include <memory>

class HomomorphismModel
{
    private:
        struct Imp;
        std::unique_ptr<Imp> _imp;

        auto _build_exact_path_graphs(std::vector<SVOBitset> & graph_rows, unsigned size, unsigned & idx,
                unsigned number_of_exact_path_graphs, bool directed) -> void;

        auto _build_distance3_graphs(std::vector<SVOBitset> & graph_rows, unsigned size, unsigned & idx) -> void;

        auto _build_k4_graphs(std::vector<SVOBitset> & graph_rows, unsigned size, unsigned & idx) -> void;

        auto _check_degree_compatibility(
                int p,
                int t,
                unsigned graphs_to_consider,
                std::vector<std::vector<std::vector<int> > > & patterns_ndss,
                std::vector<std::vector<std::optional<std::vector<int> > > > & targets_ndss,
                bool do_not_do_nds_yet
                ) const -> bool;

        auto _check_loop_compatibility(int p, int t) const -> bool;

        auto _check_label_compatibility(int p, int t) const -> bool;

        auto _bigraph_link_match(
                const std::pair<bool, std::vector<int> > & pattern_hyperedge,
                const std::pair<bool, std::vector<int> > & target_hyperedge,
                const VertexToVertexMapping & mapping) const -> bool;

        auto _backtracking_open_link_matching(
                std::vector<std::pair<bool, std::vector<int> > > & pattern_links,
                std::vector<std::pair<bool, std::vector<int> > > & target_links,
                const VertexToVertexMapping & mapping) const -> bool;

        auto _closed_link_matching(
                const VertexToVertexMapping & mapping) const -> std::vector<VertexToVertexMapping>;

        auto _check_closed_link_compatibility(int p, int t) const -> bool;

        auto _check_open_directed_hyperedges(int p, int t) const -> bool;
        
        auto _check_bigraph_equality_regions(int p, int t) const -> bool;

        auto _check_bigraph_equality_links(int p, int t) const -> bool;
        
        auto _check_bigraph_degree_compatibility(int p, int t) const -> bool;

    public:
        using PatternAdjacencyBitsType = uint8_t;

        const unsigned max_graphs;
        unsigned pattern_size, target_size;
        unsigned pattern_link_count, target_link_count;

        auto has_less_thans() const -> bool;
        std::vector<std::pair<unsigned, unsigned> > pattern_less_thans_in_convenient_order;

        HomomorphismModel(const InputGraph & target, const InputGraph & pattern, const HomomorphismParams & params);
        ~HomomorphismModel();

        auto pattern_vertex_for_proof(int v) const -> NamedVertex;
        auto target_vertex_for_proof(int v) const -> NamedVertex;

        auto prepare() -> bool;

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

        auto initialise_domains(std::vector<HomomorphismDomain> & domains) const -> bool;

        auto is_pattern_anchor(int p) const -> bool;
        auto check_extra_bigraph_constraints(const VertexToVertexMapping & mapping) const -> bool;
};

#endif
