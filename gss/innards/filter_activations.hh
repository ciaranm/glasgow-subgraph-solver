#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_FILTER_ACTIVATIONS_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_FILTER_ACTIVATIONS_HH 1

#include <list>
#include <string>
#include <vector>

namespace gss::innards
{
    // How much each filter actually removed during a solve. This exists for the option
    // sweep (gss/option_sweep_test.cc): without it the sweep cannot tell "this filter is
    // correct here" from "this filter did nothing here", and a great many of its cells are
    // vacuous by construction -- the k4 graph removes nothing on a triangle-free instance,
    // clique-size constraints remove nothing when a pattern has no clique bigger than an
    // edge, and the traits layer *silently* switches whole filters off for some option
    // combinations rather than throwing.
    //
    // Only maintained when HomomorphismParams::record_filter_activations is set, and read
    // out of HomomorphismResult::extra_stats rather than through the API. The searcher's
    // counters are behind a branch in its forward-checking loop because attributing a
    // removal to one graph pair costs a popcount per graph pair; the model's are on the
    // once-per-solve domain initialisation path, where the cost does not matter.
    struct FilterActivations
    {
        // Initial domain construction: (pattern vertex, target vertex) pairs rejected,
        // attributed to the first filter that rejected the pair, since the checks
        // short-circuit. `degree` is indexed by graph pair, and counts a rejection by the
        // degree bound on that graph (graph 0 is the original graph, the rest supplemental).
        unsigned long long vertex_labels = 0, loops = 0, nds = 0, cliques = 0;
        std::vector<unsigned long long> degree;

        // Forward checking during search: values removed from a domain, by graph pair, plus
        // the edge-label check that follows them.
        std::vector<unsigned long long> search;
        unsigned long long search_edge_labels = 0;

        explicit FilterActivations(unsigned max_graphs);

        // Did any filter that an option can turn off remove anything? Adjacency in the
        // original graph, vertex labels and loops are part of the problem rather than
        // optional filtering, so they do not count towards this.
        [[nodiscard]] auto any_optional_filter_fired() const -> bool;

        // One "key = tokens" line each, in the style of the other extra stats, or nothing
        // when that half was never filled in. Values that are per graph pair are written as
        // a comma-separated list in graph order.
        auto add_initial_stats(std::list<std::string> & stats) const -> void;
        auto add_search_stats(std::list<std::string> & stats) const -> void;
    };
}

#endif
