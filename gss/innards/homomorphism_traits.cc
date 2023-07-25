#include <gss/innards/homomorphism_traits.hh>

using namespace gss;
using namespace gss::innards;

auto gss::innards::supports_exact_path_graphs(const HomomorphismParams & params) -> bool
{
    return (! params.no_supplementals) && (params.injectivity != Injectivity::NonInjective);
}

auto gss::innards::supports_distance2_graphs(const HomomorphismParams & params) -> bool
{
    // exact path graphs are better
    return (! params.no_supplementals) && (! supports_exact_path_graphs(params));
}

auto gss::innards::supports_k4_graphs(const HomomorphismParams & params) -> bool
{
    return (! params.no_supplementals) && params.k4 && (params.injectivity != Injectivity::NonInjective);
}

auto gss::innards::supports_distance3_graphs(const HomomorphismParams & params) -> bool
{
    return (! params.no_supplementals) && params.distance3 && (params.injectivity == Injectivity::Injective);
}

auto gss::innards::might_have_watches(const HomomorphismParams & params) -> bool
{
    return params.restarts_schedule->might_restart();
}

auto gss::innards::is_nonshrinking(const HomomorphismParams & params) -> bool
{
    return params.injectivity == Injectivity::Injective;
}

auto gss::innards::degree_and_nds_are_preserved(const HomomorphismParams & params) -> bool
{
    return params.injectivity != Injectivity::NonInjective;
}

auto gss::innards::degree_and_nds_are_exact(const HomomorphismParams & params, unsigned pattern_size, unsigned target_size) -> bool
{
    return params.induced && pattern_size == target_size;
}

auto gss::innards::global_degree_is_preserved(const HomomorphismParams & params) -> bool
{
    return params.injectivity == Injectivity::Injective;
}

auto gss::innards::can_use_clique(const HomomorphismParams & params) -> bool
{
    return (! params.count_solutions) && (! params.lackey) && params.clique_detection && (! params.proof_options);
}
