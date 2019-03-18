/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_TRAITS_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_TRAITS_HH 1

#include "homomorphism.hh"

auto can_strip_isolated_vertices(const HomomorphismParams & params) -> bool
{
    return (! params.induced) && (! params.count_solutions) && params.remove_isolated_vertices;
}

auto supports_exact_path_graphs(const HomomorphismParams & params) -> bool
{
    return params.injectivity == Injectivity::Injective;
}

auto might_have_watches(const HomomorphismParams & params) -> bool
{
    return ! params.count_solutions;
}

auto is_nonshrinking(const HomomorphismParams & params) -> bool
{
    return params.injectivity == Injectivity::Injective;
}

auto degree_and_nds_are_preserved(const HomomorphismParams & params) -> bool
{
    return params.injectivity != Injectivity::NonInjective;
}

auto can_use_clique(const HomomorphismParams & params) -> bool
{
    return (! params.count_solutions) && params.clique_detection;
}

#endif
