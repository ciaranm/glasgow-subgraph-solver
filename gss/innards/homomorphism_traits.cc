#include <gss/innards/homomorphism_traits.hh>

using namespace gss;
using namespace gss::innards;

// Local injectivity together with a self-loop makes the loop-stripped, injective-style
// supplemental and degree/NDS reasoning unsound (issue #58): a neighbour can map onto a
// target self-loop, which the stripped rows and the "distinct images" counting miss. In
// that one case fall back to adjacency + local-injectivity propagation, which is correct.
static auto loop_breaks_filtering(const gss::HomomorphismParams & params, bool has_loops) -> bool
{
    return params.injectivity == gss::Injectivity::LocallyInjective && has_loops;
}

auto gss::innards::supports_exact_path_graphs(const HomomorphismParams & params, bool has_loops) -> bool
{
    return (! params.no_supplementals) && (params.injectivity != Injectivity::NonInjective) && (! loop_breaks_filtering(params, has_loops));
}

auto gss::innards::supports_distance2_graphs(const HomomorphismParams & params, bool has_loops) -> bool
{
    // exact path graphs are better
    return (! params.no_supplementals) && (! supports_exact_path_graphs(params, has_loops)) && (! loop_breaks_filtering(params, has_loops));
}

auto gss::innards::supports_k4_graphs(const HomomorphismParams & params, bool has_loops) -> bool
{
    return (! params.no_supplementals) && params.k4 && (params.injectivity != Injectivity::NonInjective) && (! loop_breaks_filtering(params, has_loops));
}

auto gss::innards::supports_distance3_graphs(const HomomorphismParams & params) -> bool
{
    return (! params.no_supplementals) && params.distance3 && (params.injectivity == Injectivity::Injective);
}

auto gss::innards::might_have_watches(const HomomorphismParams & params) -> bool
{
    // Watches (the nogood store) are only worth maintaining if the search can restart --
    // otherwise nothing ever consults a posted nogood. Staging always performs at least one
    // restart (the Stage-1 -> Stage-2 transition), independently of the user's schedule, so
    // it must enable watches even under a no-restart schedule: the transition's
    // restart-resumption nogoods are what stop Stage 2 re-exploring (and, when counting,
    // re-counting) the region Stage 1 already searched.
    return params.staged || params.restarts_schedule->might_restart();
}

auto gss::innards::is_nonshrinking(const HomomorphismParams & params) -> bool
{
    return params.injectivity == Injectivity::Injective;
}

auto gss::innards::degree_and_nds_are_preserved(const HomomorphismParams & params, bool has_loops) -> bool
{
    return (params.injectivity != Injectivity::NonInjective) && (! loop_breaks_filtering(params, has_loops));
}

auto gss::innards::degree_and_nds_are_exact(const HomomorphismParams & params, unsigned pattern_size, unsigned target_size) -> bool
{
    // Degrees (and neighbourhood degree sequences) must match *exactly* only when the
    // mapping is forced to be a bijection -- i.e. an induced graph isomorphism. That needs
    // full injectivity: an injection between equal-sized vertex sets is onto. Under merely
    // local injectivity the mapping can collide, so a degree-d pattern vertex may map to a
    // higher-degree target (only deg(image) >= d is required), and demanding equality
    // wrongly prunes those (issue #58).
    return params.induced && pattern_size == target_size && params.injectivity == Injectivity::Injective;
}

auto gss::innards::global_degree_is_preserved(const HomomorphismParams & params) -> bool
{
    return params.injectivity == Injectivity::Injective;
}

auto gss::innards::can_use_clique(const HomomorphismParams & params) -> bool
{
    return (! params.count_solutions) && params.clique_detection && (! params.proof_options);
}
