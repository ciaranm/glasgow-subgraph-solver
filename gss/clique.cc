#include <gss/clique.hh>
#include <gss/configuration.hh>
#include <gss/innards/proof.hh>
#include <gss/innards/svo_bitset.hh>
#include <gss/innards/watches.hh>

#include <algorithm>
#include <list>
#include <memory>
#include <numeric>
#include <random>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

using namespace gss;
using namespace gss::innards;

using std::conditional_t;
using std::find;
using std::iota;
using std::is_same;
using std::list;
using std::make_shared;
using std::make_tuple;
using std::max_element;
using std::move;
using std::mt19937;
using std::pair;
using std::reverse;
using std::shared_ptr;
using std::sort;
using std::string_view;
using std::swap;
using std::to_string;
using std::vector;

namespace
{
    enum class SearchResult
    {
        Aborted,
        Restart,
        Complete,
        DecidedTrue
    };

    struct Incumbent
    {
        unsigned value = 0;
        vector<int> c;

        auto update(const vector<int> & new_c, unsigned long long & find_nodes, unsigned long long & prove_nodes) -> void
        {
            if (new_c.size() > value) {
                find_nodes += prove_nodes;
                prove_nodes = 0;
                value = new_c.size();
                c = new_c;
            }
        }
    };

    /**
     * What proof logging needs from a filter pass, kept per search node because the
     * bound that uses it is only logged once the node's branches are done with. Built
     * only when we are actually writing a proof.
     */
    struct FilterProofRecord
    {
        vector<vector<int>> unconflicted_classes;
        vector<CliqueConflict> conflicts;
    };

    template <typename EntryType_>
    struct FlatWatchTable
    {
        vector<EntryType_> data;

        EntryType_ & operator[](int x)
        {
            return data[x];
        }
    };

    struct CliqueRunner
    {
        const CliqueParams & params;
        Incumbent incumbent;
        shared_ptr<Proof> proof;

        int size;
        vector<SVOBitset> adj, connected_table;
        vector<int> order, invorder;

        Watches<int, FlatWatchTable> watches;

        mt19937 global_rand;

        std::unique_ptr<int[]> space;

        // Scratch for the branching-set filter. The colour classes below the pruning
        // threshold are held as intrusive singly linked lists, so that moving a vertex
        // between classes is O(1) and walking a class can bail out early. Only live
        // during a single filter pass, which finishes before the node recurses, so one
        // copy for the whole search is enough and nothing is allocated per node.
        vector<int> cl_next, cl_head;
        vector<unsigned char> cl_forbidden;

        CliqueRunner(const InputGraph & g, const CliqueParams & p) :
            params(p),
            size(g.size()),
            adj(g.size(), SVOBitset{unsigned(size), 0}),
            order(size),
            invorder(size)
        {
            if (p.proof_options)
                proof = make_shared<Proof>(*p.proof_options);
            else if (p.extend_proof)
                proof = p.extend_proof;

            if (proof && ! proof->has_clique_model() && ! params.proof_is_for_hom) {
                for (int q = 0; q < g.size(); ++q)
                    proof->create_binary_variable(q, [&](int v) { return g.vertex_name(v); });

                proof->create_objective(g.size(), params.decide);

                for (int p = 0; p < g.size(); ++p)
                    for (int q = 0; q < p; ++q)
                        if (! g.adjacent(p, q))
                            proof->create_non_edge_constraint(pair{p, g.vertex_name(p)}, pair{q, g.vertex_name(q)});

                proof->finalise_model();
            }

            if (params.restarts_schedule->might_restart())
                watches.table.data.resize(g.size());

            // populate our order with every vertex initially
            iota(order.begin(), order.end(), 0);

            // pre-calculate degrees
            vector<int> degrees;
            degrees.resize(size);
            g.for_each_edge([&](int f, int t, string_view) { if (f != t) ++degrees[f]; });

            // Workspace for the search: two int[size] arrays per recursion level. The
            // search never recurses deeper than the largest clique, which is at most
            // one more than the largest degree, so size the workspace to that rather
            // than to the number of vertices (which would be quadratic in memory, and
            // overflows an int product beyond ~32k vertices -- issue #39). size_t
            // arithmetic keeps it correct for dense graphs where the degree is large.
            int max_degree = degrees.empty() ? 0 : *max_element(degrees.begin(), degrees.end());
            space = std::make_unique<int[]>(std::size_t(2) * size * (max_degree + 2));

            if (params.filter != CliqueFilter::None && ! params.connected) {
                cl_next.resize(size);
                cl_head.resize(size + 1);
                cl_forbidden.resize(size + 1);
            }

            // sort on degree
            if (! params.input_order)
                sort(order.begin(), order.end(),
                    [&](int a, int b) { return (degrees[a] > degrees[b] || (degrees[a] == degrees[b] && a < b)); });

            for (unsigned i = 0; i < order.size(); ++i)
                invorder[order[i]] = i;

            // Loops are irrelevant to cliques (a vertex is never its own clique
            // neighbour). Including them would leave a vertex's own bit set in its
            // adjacency, letting it be re-selected during search and recursing past
            // the workspace bound — a crash on looped inputs (issue #38).
            g.for_each_edge([&](int f, int t, string_view) { if (f != t) adj[invorder[f]].set(invorder[t]); });

            if (params.connected) {
                connected_table.resize(size);
                for (int v = 0; v < size; ++v)
                    connected_table[v] = params.connected(order.at(v), [&](int x) { return invorder.at(x); });
            }
        }

        auto colour_class_order(
            const SVOBitset & p,
            int * p_order,
            int * p_bounds,
            int & p_end) -> void
        {
            SVOBitset p_left = p; // not coloured yet
            unsigned colour = 0;  // current colour
            p_end = 0;

            // while we've things left to colour
            while (p_left.any()) {
                // next colour
                ++colour;
                // things that can still be given this colour
                SVOBitset q = p_left;

                // while we can still give something this colour
                while (q.any()) {
                    // first thing we can colour
                    int v = q.find_first();
                    p_left.reset(v);
                    q.reset(v);

                    // can't give anything adjacent to this the same colour
                    q.intersect_with_complement(adj[v]);

                    // record in result
                    p_bounds[p_end] = colour;
                    p_order[p_end] = v;
                    ++p_end;
                }
            }
        }

        auto connected_colour_class_order(
            const SVOBitset & p,
            const SVOBitset & a,
            int * p_order,
            int * p_bounds,
            int & p_end) -> void
        {
            unsigned colour = 0; // current colour
            p_end = 0;

            SVOBitset p_left = p; // not coloured yet
            p_left.intersect_with_complement(a);

            // while we've things left to colour
            while (p_left.any()) {
                // next colour
                ++colour;
                // things that can still be given this colour
                SVOBitset q = p_left;

                // while we can still give something this colour
                while (q.any()) {
                    // first thing we can colour
                    int v = q.find_first();
                    p_left.reset(v);
                    q.reset(v);

                    // can't give anything adjacent to this the same colour
                    q.intersect_with_complement(adj[v]);

                    // record in result
                    p_bounds[p_end] = colour;
                    p_order[p_end] = v;
                    ++p_end;
                }
            }

            p_left = p;
            p_left &= a;

            // while we've things left to colour
            while (p_left.any()) {
                // next colour
                ++colour;
                // things that can still be given this colour
                SVOBitset q = p_left;

                // while we can still give something this colour
                while (q.any()) {
                    // first thing we can colour
                    int v = q.find_first();
                    p_left.reset(v);
                    q.reset(v);

                    // can't give anything adjacent to this the same colour
                    q.intersect_with_complement(adj[v]);

                    // record in result
                    p_bounds[p_end] = colour;
                    p_order[p_end] = v;
                    ++p_end;
                }
            }
        }

        auto colour_class_order_2df(
            const SVOBitset & p,
            int * p_order,
            int * p_bounds,
            int * defer,
            int & p_end) -> void
        {
            SVOBitset p_left = p; // not coloured yet
            unsigned colour = 0;  // current colour
            p_end = 0;

            unsigned d = 0; // number deferred

            // while we've things left to colour
            while (p_left.any()) {
                // next colour
                ++colour;
                // things that can still be given this colour
                SVOBitset q = p_left;

                // while we can still give something this colour
                unsigned number_with_this_colour = 0;
                while (q.any()) {
                    // first thing we can colour
                    int v = q.find_first();
                    p_left.reset(v);
                    q.reset(v);

                    // can't give anything adjacent to this the same colour
                    q.intersect_with_complement(adj[v]);

                    // record in result
                    p_bounds[p_end] = colour;
                    p_order[p_end] = v;
                    ++p_end;
                    ++number_with_this_colour;
                }

                if (1 == number_with_this_colour) {
                    --p_end;
                    --colour;
                    defer[d++] = p_order[p_end];
                }
            }

            // handle deferred singletons
            for (unsigned n = 0; n < d; ++n) {
                ++colour;
                p_order[p_end] = defer[n];
                p_bounds[p_end] = colour;
                ++p_end;
            }
        }

        auto colour_class_order_sorted(
            const SVOBitset & p,
            int * p_order,
            int * p_bounds,
            int & p_end) -> void
        {
            SVOBitset p_left = p; // not coloured yet
            unsigned colour = 0;  // current colour
            p_end = 0;

            vector<int> p_order_prelim(size);
            vector<int> colour_sizes(size);
            vector<int> colour_start(size);
            vector<int> sorted_order(size);

            // while we've things left to colour
            while (p_left.any()) {
                colour_start[colour] = p_end;
                colour_sizes[colour] = 0;

                // next colour
                ++colour;
                // things that can still be given this colour
                SVOBitset q = p_left;

                // while we can still give something this colour
                while (q.any()) {
                    // first thing we can colour
                    int v = q.find_first();
                    p_left.reset(v);
                    q.reset(v);

                    // can't give anything adjacent to this the same colour
                    q.intersect_with_complement(adj[v]);

                    // record in result
                    p_order_prelim[p_end] = v;
                    ++p_end;
                    ++colour_sizes[colour - 1];
                }
            }

            // sort
            iota(sorted_order.begin(), sorted_order.begin() + colour, 0);
            sort(sorted_order.begin(), sorted_order.begin() + colour, [&](int a, int b) {
                return make_tuple(colour_sizes[b], a) < make_tuple(colour_sizes[a], b);
            });

            // copy out
            int p_end2 = 0;
            for (unsigned c = 0; c < colour; ++c) {
                for (int v = colour_start[sorted_order[c]]; v < colour_start[sorted_order[c]] + colour_sizes[sorted_order[c]]; ++v) {
                    p_bounds[p_end2] = c + 1;
                    p_order[p_end2] = p_order_prelim[v];
                    ++p_end2;
                }
            }
        }

        /**
         * Try to show that vertex v need not be branched on, given the colour classes
         * below the pruning threshold.
         *
         * Both of the tests here hang off the same quantity, C_k1 & N(v), which is why
         * they are fused: San Segundo et al's FILTER_RECOL_INFRACHROM. If v has no
         * neighbour in some class, it simply belongs there (one-move recolouring). If it
         * has exactly one neighbour w there, then either w can be moved out to another
         * class and v can take its place (double-move recolouring, Tomita's Re-NUMBER),
         * or no vertex of a third class is adjacent to both v and w -- in which case
         * ({v}, C_k1, C_k2) admits no triangle with one vertex in each, and contributes
         * two rather than three to the bound.
         *
         * The classes are walked rather than intersected as bitsets on purpose. The
         * common case is failure, and failure exits after finding a second neighbour of
         * v, which on a dense graph happens after a couple of tests; a bitset
         * intersection would pay for every word every time. Prosser's observation that a
         * vertex of colour k has a neighbour in every class below k means the first
         * neighbour is always there to be found.
         */
        auto filter_one(int v, unsigned k_min, FilterProofRecord * rec) -> bool
        {
            for (unsigned k1 = 1; k1 < k_min; ++k1) {
                if (cl_forbidden[k1])
                    continue;

                // how many neighbours does v have in class k1? we only care whether it
                // is none, exactly one, or more than one, so stop counting at two
                int w = -1, w_prev = -1, prev = -1, hits = 0;
                for (int u = cl_head[k1]; u != -1; prev = u, u = cl_next[u])
                    if (adj[v].test(u)) {
                        if (++hits > 1)
                            break;
                        w = u;
                        w_prev = prev;
                    }

                if (hits > 1)
                    continue;

                if (0 == hits) {
                    // one-move recolouring: nothing in this class conflicts with v, so
                    // v belongs in it and is bounded by k1 < k_min
                    cl_next[v] = cl_head[k1];
                    cl_head[k1] = v;
                    return true;
                }

                for (unsigned k2 = 1; k2 < k_min; ++k2) {
                    if (k2 == k1 || cl_forbidden[k2])
                        continue;

                    // one walk answers both questions: does w have a neighbour in class
                    // k2 at all, and does it have one that is also adjacent to v?
                    bool w_has_neighbour = false, common_neighbour = false;
                    for (int u = cl_head[k2]; u != -1; u = cl_next[u])
                        if (adj[w].test(u)) {
                            w_has_neighbour = true;
                            if (adj[v].test(u)) {
                                common_neighbour = true;
                                break;
                            }
                        }

                    if (common_neighbour)
                        continue;

                    if (! w_has_neighbour) {
                        // double-move recolouring: w is free to join class k2, which
                        // leaves class k1 open for v
                        if (-1 == w_prev)
                            cl_head[k1] = cl_next[w];
                        else
                            cl_next[w_prev] = cl_next[w];
                        cl_next[w] = cl_head[k2];
                        cl_head[k2] = w;
                        cl_next[v] = cl_head[k1];
                        cl_head[k1] = v;
                        return true;
                    }

                    if (CliqueFilter::InfraChromatic != params.filter)
                        continue;

                    // ({v}, C_k1, C_k2) is conflicting. Both colours are now spent:
                    // the bound is only sound if the conflicts we collect are disjoint,
                    // because it is really one at-most-two constraint over the union of
                    // the three classes replacing three at-most-ones.
                    if (rec) {
                        CliqueConflict cf;
                        cf.filtered_vertex = order[v];
                        for (auto k : {k1, k2})
                            for (int u = cl_head[k]; u != -1; u = cl_next[u]) {
                                (k == k1 ? cf.class1 : cf.class2).push_back(order[u]);
                                // the part of the two classes v can see is an independent
                                // set: it is just w from k1, and nothing in k2 adjacent to
                                // v is adjacent to w
                                (adj[v].test(u) ? cf.independent_set : cf.non_neighbours).push_back(order[u]);
                            }
                        rec->conflicts.push_back(move(cf));
                    }

                    cl_forbidden[k1] = 1;
                    cl_forbidden[k2] = 1;
                    return true;
                }
            }

            return false;
        }

        /**
         * Thin out the branching set. Vertices whose colour is below k_min are already
         * bounded; for each of the rest, in increasing colour order, see whether
         * filter_one can dispose of it.
         *
         * Once this is done, everything still in p when the bound finally fires is
         * covered by the (possibly recoloured) classes below k_min together with one
         * singleton per infra-chromatic conflict, and each conflict knocks one off the
         * bound. That is why the colours a conflict uses are marked forbidden: the sum
         * is only k_min - 1 if the conflicts are over disjoint sets of classes.
         *
         * Returns the lowest position filtered, or p_end if none were.
         */
        auto recolour_and_filter(
            const int * p_order,
            int * p_bounds,
            int p_end,
            unsigned k_min,
            FilterProofRecord * rec) -> int
        {
            for (unsigned k = 1; k < k_min; ++k) {
                cl_head[k] = -1;
                cl_forbidden[k] = 0;
            }

            int first_branch = 0;
            while (first_branch < p_end && unsigned(p_bounds[first_branch]) < k_min)
                ++first_branch;

            // build the class lists, in increasing position order so that the search is
            // deterministic and matches the colouring
            for (int i = first_branch - 1; i >= 0; --i) {
                int v = p_order[i];
                cl_next[v] = cl_head[p_bounds[i]];
                cl_head[p_bounds[i]] = v;
            }

            int min_filtered = p_end;
            for (int i = first_branch; i < p_end; ++i)
                if (filter_one(p_order[i], k_min, rec)) {
                    // negating the colour marks the position as not worth branching on.
                    // The vertex deliberately stays in p: all we have shown is that we
                    // need not branch on it here, not that no clique uses it, and a
                    // clique through some other branching vertex may still need it.
                    p_bounds[i] = -p_bounds[i];
                    if (p_end == min_filtered)
                        min_filtered = i;
                }

            // the classes as they now stand, minus the ones spent on a conflict: those
            // reach the proof through the conflicts instead, at two apiece rather than one
            if (rec)
                for (unsigned k = 1; k < k_min; ++k)
                    if (! cl_forbidden[k] && -1 != cl_head[k]) {
                        rec->unconflicted_classes.emplace_back();
                        for (int u = cl_head[k]; u != -1; u = cl_next[u])
                            rec->unconflicted_classes.back().push_back(order[u]);
                    }

            return min_filtered;
        }

        auto post_nogood(
            const vector<int> & c)
        {
            Nogood<int> nogood;
            nogood.literals.assign(c.begin(), c.end());
            watches.post_nogood(move(nogood));
        }

        auto unpermute(
            const vector<int> & v) -> vector<int>
        {
            vector<int> result;
            for (auto & w : v)
                result.push_back(order[w]);
            return result;
        }

        auto unpermute_and_finish(
            vector<int> & v) -> vector<pair<int, bool>>
        {
            vector<pair<int, bool>> result;
            for (auto & w : v)
                result.emplace_back(order[w], true);
            for (int w = 0; w < size; ++w)
                if (result.end() == find_if(result.begin(), result.end(), [&](auto & x) { return x.first == w; }))
                    result.emplace_back(w, false);
            return result;
        }

        template <bool connected_>
        auto expand(
            int depth,
            unsigned long long & nodes,
            unsigned long long & find_nodes,
            unsigned long long & prove_nodes,
            vector<int> & c,
            SVOBitset & p,
            conditional_t<connected_, const SVOBitset &, int> a,
            int spacepos) -> SearchResult
        {
            ++nodes;
            ++prove_nodes;

            // initial colouring
            int * p_order = &space[spacepos];
            int * p_bounds = &space[spacepos + size];

            int p_end = 0;

            if constexpr (connected_) {
                if (! c.empty())
                    connected_colour_class_order(p, a, p_order, p_bounds, p_end);
                else
                    colour_class_order(p, p_order, p_bounds, p_end);
            }
            else {
                switch (params.colour_class_order) {
                case ColourClassOrder::ColourOrder: colour_class_order(p, p_order, p_bounds, p_end); break;
                case ColourClassOrder::SingletonsFirst: colour_class_order_2df(p, p_order, p_bounds, &space[spacepos + 2 * size], p_end); break;
                case ColourClassOrder::Sorted: colour_class_order_sorted(p, p_order, p_bounds, p_end); break;
                }
            }

            // Anything coloured below k_min is already bounded; everything at or above it
            // is a vertex we would otherwise have to branch on. See how many of those we
            // can dispose of without branching. Filtered positions come back with their
            // colour negated.
            int min_filtered_pos = p_end;
            unsigned filter_k_min = 0;
            std::unique_ptr<FilterProofRecord> filter_record;
            if constexpr (! connected_) {
                if (CliqueFilter::None != params.filter) {
                    unsigned k_min = incumbent.value >= c.size() ? incumbent.value - c.size() + 1 : 1;
                    if (k_min >= 2) {
                        filter_k_min = k_min;
                        if (proof)
                            filter_record = std::make_unique<FilterProofRecord>();
                        min_filtered_pos = recolour_and_filter(p_order, p_bounds, p_end, k_min, filter_record.get());
                    }
                }
            }

            // for each v in p... (v comes later)
            for (int n = p_end - 1; n >= 0; --n) {
                // bound, timeout or early exit?
                if (params.timeout->should_abort())
                    return SearchResult::Aborted;

                // filtered out by recolouring or an infra-chromatic conflict: do not
                // branch on it, but leave it in p for the branches that come after
                if (p_bounds[n] < 0)
                    continue;

                if (c.size() + p_bounds[n] <= incumbent.value) {
                    if (proof && filter_record) {
                        // The classes below the filter threshold are the recoloured ones,
                        // and the ones spent on a conflict are left out because the
                        // conflicts cover them. Above the threshold the colouring is
                        // untouched, except that filtered vertices are skipped: they have
                        // already been accounted for, either by having been recoloured
                        // into a class below the threshold or by being a conflict's
                        // singleton.
                        auto colour_classes = filter_record->unconflicted_classes;
                        int previous_colour = 0;
                        for (int v = 0; v <= n; ++v) {
                            if (p_bounds[v] < 0 || unsigned(p_bounds[v]) < filter_k_min)
                                continue;
                            if (p_bounds[v] != previous_colour) {
                                colour_classes.emplace_back();
                                previous_colour = p_bounds[v];
                            }
                            colour_classes.back().push_back(order[p_order[v]]);
                        }
                        proof->colour_bound(colour_classes, filter_record->conflicts);
                    }
                    else if (proof) {
                        vector<vector<int>> colour_classes;
                        for (int v = 0; v <= n; ++v) {
                            if (0 == v || p_bounds[v - 1] != p_bounds[v])
                                colour_classes.emplace_back();
                            colour_classes.back().push_back(order[p_order[v]]);
                        }
                        proof->colour_bound(colour_classes);
                    }
                    break;
                }

                // if we've used k colours to colour k vertices, it's a clique. this isn't (I think?) a
                // valid shortcut in the connected case.
                if constexpr (! connected_) {
                    if (n < min_filtered_pos && p_bounds[n] == n + 1) {
                        auto c_save = c;
                        for (; n >= 0; --n)
                            c.push_back(p_order[n]);
                        incumbent.update(c, find_nodes, prove_nodes);

                        if (proof && ! params.decide) {
                            proof->start_level(0);
                            proof->new_incumbent(unpermute_and_finish(c));
                            proof->start_level(depth + 1);
                        }

                        if ((params.decide && incumbent.value >= *params.decide) ||
                            (params.stop_after_finding && incumbent.value >= *params.stop_after_finding)) {
                            if (proof)
                                proof->post_solution(unpermute(c));

                            return SearchResult::DecidedTrue;
                        }

                        c = move(c_save);

                        break;
                    }
                }

                auto v = p_order[n];

                if constexpr (connected_) {
                    if ((! c.empty()) && (! a.test(v))) {
                        // none of the remaining vertices can give a connected underlying graph
                        if (proof) {
                            auto c_unpermuted = unpermute(c);
                            for (int v = 0; v <= n; ++v)
                                proof->not_connected_in_underlying_graph(unpermute(c), order[p_order[v]]);

                            proof->start_level(depth);
                            proof->backtrack_from_binary_variables(unpermute(c));
                            proof->forget_level(depth + 1);
                        }

                        break;
                    }
                }

                // consider taking v
                c.push_back(v);

                if (params.decide || params.stop_after_finding) {
                    if ((params.decide && incumbent.value >= *params.decide) ||
                        (params.stop_after_finding && incumbent.value >= *params.stop_after_finding)) {
                        if (proof)
                            proof->post_solution(unpermute(c));

                        return SearchResult::DecidedTrue;
                    }
                }
                else {
                    if (proof && c.size() > incumbent.value && ! params.proof_is_for_hom) {
                        proof->start_level(0);
                        proof->new_incumbent(unpermute_and_finish(c));
                        proof->start_level(depth + 1);
                    }
                    incumbent.update(c, find_nodes, prove_nodes);
                }

                // filter p to contain vertices adjacent to v
                SVOBitset new_p = p;
                new_p &= adj[v];

                if (params.restarts_schedule->might_restart())
                    watches.propagate(
                        v,
                        [&](int literal) { return c.end() == find(c.begin(), c.end(), literal); },
                        [&](int literal) { new_p.reset(literal); });

                if (proof)
                    proof->start_level(depth + 1);

                if (new_p.any()) {
                    auto new_a = a;

                    if constexpr (connected_) {
                        new_a |= connected_table[v];
                    }

                    switch (expand<connected_>(depth + 1, nodes, find_nodes, prove_nodes, c, new_p, new_a, spacepos + 2 * size)) {
                    case SearchResult::Aborted:
                        return SearchResult::Aborted;

                    case SearchResult::DecidedTrue:
                        return SearchResult::DecidedTrue;

                    case SearchResult::Complete:
                        break;

                    case SearchResult::Restart:
                        // restore assignments before posting nogoods, it's easier
                        c.pop_back();

                        // post nogoods for everything we've done so far. A filtered
                        // vertex is not "done": all we know is that the classes below
                        // k_min plus the filtered vertices cannot beat the incumbent
                        // between them, which says nothing until the rest of the
                        // branching set has been dealt with, and here it has not.
                        for (int m = p_end - 1; m > n; --m) {
                            if (p_bounds[m] < 0)
                                continue;
                            c.push_back(p_order[m]);
                            post_nogood(c);
                            c.pop_back();
                        }

                        return SearchResult::Restart;
                    }
                }

                if (proof) {
                    proof->start_level(depth);
                    proof->backtrack_from_binary_variables(unpermute(c));
                    proof->forget_level(depth + 1);
                }

                // now consider not taking v
                c.pop_back();
                p.reset(v);
            }

            params.restarts_schedule->did_a_backtrack();
            if (params.restarts_schedule->should_restart()) {
                post_nogood(c);
                return SearchResult::Restart;
            }
            else
                return SearchResult::Complete;
        }

        template <bool connected_>
        auto run() -> CliqueResult
        {
            CliqueResult result;

            if (params.decide)
                incumbent.value = *params.decide - 1;

            // do the search
            bool done = false;
            unsigned number_of_restarts = 0;

            SVOBitset p{unsigned(size), 0};
            for (int i = 0; i < size; ++i)
                p.set(i);

            while (! done) {
                ++number_of_restarts;

                // start watching new nogoods
                done = watches.apply_new_nogoods(
                    [&](int literal) { p.reset(literal); });

                if (done)
                    break;

                watches.clear_new_nogoods();

                auto new_p = p;
                vector<int> c;
                conditional_t<connected_, SVOBitset, int> a{};
                if constexpr (connected_)
                    a = SVOBitset{unsigned(size), 0};

                switch (expand<connected_>(params.proof_is_for_hom ? 1 : 0, result.nodes, result.find_nodes, result.prove_nodes, c, new_p, a, 0)) {
                case SearchResult::Complete:
                    done = true;
                    break;

                case SearchResult::DecidedTrue:
                    done = true;
                    break;

                case SearchResult::Aborted:
                    done = true;
                    break;

                case SearchResult::Restart:
                    break;
                }

                params.restarts_schedule->did_a_restart();
            }

            if (params.restarts_schedule->might_restart())
                result.extra_stats.emplace_back("restarts = " + to_string(number_of_restarts));

            if (proof && params.decide && incumbent.c.empty() && ! params.proof_is_for_hom)
                proof->finish_unsat_proof();
            else if (proof && ! params.decide && params.adjust_objective_for_mcs)
                proof->finish_optimisation_proof(*params.adjust_objective_for_mcs - incumbent.c.size());
            else if (proof && ! params.decide && ! params.proof_is_for_hom)
                proof->finish_optimisation_proof(size - incumbent.c.size());

            result.clique.clear();
            for (auto & v : incumbent.c)
                result.clique.insert(order[v]);

            return result;
        }
    };
}

auto gss::solve_clique_problem(const InputGraph & graph, const CliqueParams & params) -> CliqueResult
{
    CliqueRunner runner{graph, params};
    return params.connected ? runner.run<true>() : runner.run<false>();
}
