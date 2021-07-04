/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "homomorphism_traits.hh"

auto supports_exact_path_graphs(const HomomorphismParams & params) -> bool
{
    return (! params.no_supplementals) && (params.injectivity != Injectivity::NonInjective);
}

auto supports_distance2_graphs(const HomomorphismParams & params) -> bool
{
    // exact path graphs are better
    return (! params.no_supplementals) && (! supports_exact_path_graphs(params));
}

auto supports_k4_graphs(const HomomorphismParams & params) -> bool
{
    return (! params.no_supplementals) && params.k4 && (params.injectivity != Injectivity::NonInjective);
}

auto supports_distance3_graphs(const HomomorphismParams & params) -> bool
{
    return (! params.no_supplementals) && params.distance3 && (params.injectivity == Injectivity::Injective);
}

auto might_have_watches(const HomomorphismParams & params) -> bool
{
    return params.restarts_schedule->might_restart();
}

auto is_nonshrinking(const HomomorphismParams & params) -> bool
{
    return params.injectivity == Injectivity::Injective;
}

auto degree_and_nds_are_preserved(const HomomorphismParams & params) -> bool
{
    return params.injectivity != Injectivity::NonInjective;
}

auto degree_and_nds_are_exact(const HomomorphismParams & params, unsigned pattern_size, unsigned target_size) -> bool
{
    return params.induced && pattern_size == target_size;
}

auto global_degree_is_preserved(const HomomorphismParams & params) -> bool
{
    return params.injectivity == Injectivity::Injective;
}

auto can_use_clique(const HomomorphismParams & params) -> bool
{
    return (! params.count_solutions) && (! params.lackey) && params.clique_detection && (! params.proof);
}

