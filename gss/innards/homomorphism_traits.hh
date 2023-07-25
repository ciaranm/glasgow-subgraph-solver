#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_TRAITS_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_TRAITS_HH 1

#include <gss/homomorphism.hh>

namespace gss::innards
{
    auto supports_exact_path_graphs(const HomomorphismParams & params) -> bool;

    auto supports_distance2_graphs(const HomomorphismParams & params) -> bool;

    auto supports_k4_graphs(const HomomorphismParams & params) -> bool;

    auto supports_distance3_graphs(const HomomorphismParams & params) -> bool;

    auto might_have_watches(const HomomorphismParams & params) -> bool;

    auto is_nonshrinking(const HomomorphismParams & params) -> bool;

    auto degree_and_nds_are_preserved(const HomomorphismParams & params) -> bool;

    auto degree_and_nds_are_exact(const HomomorphismParams & params, unsigned pattern_size, unsigned target_size) -> bool;

    auto global_degree_is_preserved(const HomomorphismParams & params) -> bool;

    auto can_use_clique(const HomomorphismParams & params) -> bool;
}

#endif
