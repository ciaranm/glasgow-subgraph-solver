#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_TRAITS_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_TRAITS_HH 1

#include <gss/homomorphism.hh>

namespace gss::innards
{
    // has_loops is whether the pattern or target has any self-loop. The supplemental-graph
    // and degree/NDS filters reason with loop-stripped rows and an injective counting
    // argument; under *local* injectivity with a loop present that argument is unsound (a
    // neighbour may map onto a target self-loop), so these are disabled for that case and
    // the search falls back to adjacency + local-injectivity propagation (issue #58).
    auto supports_exact_path_graphs(const HomomorphismParams & params, bool has_loops) -> bool;

    auto supports_distance2_graphs(const HomomorphismParams & params, bool has_loops) -> bool;

    auto supports_k4_graphs(const HomomorphismParams & params, bool has_loops) -> bool;

    auto supports_distance3_graphs(const HomomorphismParams & params) -> bool;

    auto might_have_watches(const HomomorphismParams & params) -> bool;

    auto is_nonshrinking(const HomomorphismParams & params) -> bool;

    auto degree_and_nds_are_preserved(const HomomorphismParams & params, bool has_loops) -> bool;

    auto degree_and_nds_are_exact(const HomomorphismParams & params, unsigned pattern_size, unsigned target_size) -> bool;

    auto global_degree_is_preserved(const HomomorphismParams & params) -> bool;

    // The clique-size rule (--cliques): a pattern vertex whose largest clique has k vertices
    // can only be mapped to a target vertex whose largest clique has at least k. That needs the
    // k pattern clique vertices to reach k *distinct* targets. Full injectivity gives that; so,
    // on the original graph pair, does a loopless target, because a homomorphism cannot collapse
    // two adjacent vertices onto one image without that image carrying a self-loop. Add a target
    // loop and the second argument goes, and the rule prunes solutions away (issue #91). This is
    // deliberately a different condition from the degree rule above -- same shape of counting
    // argument, different premise -- which is why it gets its own predicate.
    auto supports_clique_size_constraints(const HomomorphismParams & params, bool has_loops) -> bool;

    // The same rule on the supplemental graph pairs (--cliques-on-supplementals), where the
    // loopless argument does not apply at all: adjacency in a distance-2 graph means "within
    // distance two", and two such pattern vertices may legitimately share an image in a loopless
    // target -- build_exact_path_graphs sets that graph's diagonal precisely so that propagation
    // permits it. So this one needs injectivity, although local injectivity is enough: on the
    // supplemental graphs that actually get built, adjacency implies a common neighbour, which
    // local injectivity forces apart (issue #91).
    auto supports_clique_size_constraints_on_supplementals(const HomomorphismParams & params, bool has_loops) -> bool;

    // The clique shortcut: replace a clique pattern with a search for k pairwise-adjacent
    // target vertices. Those are distinct, so this is the same counting argument as
    // supports_clique_size_constraints() above and it needs the same premise -- full
    // injectivity, or a loopless target, without which two adjacent pattern vertices may
    // collapse onto one looped image and the reduction misses the solution (issue #94).
    auto can_use_clique(const HomomorphismParams & params, bool target_has_loops) -> bool;
}

#endif
