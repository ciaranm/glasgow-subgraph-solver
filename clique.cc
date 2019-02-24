/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "clique.hh"
#include "template_voodoo.hh"
#include "watches.hh"

#include <algorithm>
#include <list>
#include <numeric>
#include <random>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

using std::find;
using std::iota;
using std::is_same;
using std::list;
using std::make_tuple;
using std::mt19937;
using std::move;
using std::pair;
using std::reverse;
using std::sort;
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

    template <typename EntryType_>
    struct FlatWatchTable
    {
        vector<EntryType_> data;

        EntryType_ & operator[] (int x)
        {
            return data[x];
        }
    };

    template <typename BitSetType_, typename ArrayType_>
    struct CliqueRunner
    {
        const CliqueParams & params;
        Incumbent incumbent;

        int size;
        vector<BitSetType_> adj;
        vector<int> order, invorder;

        Watches<int, FlatWatchTable> watches;

        mt19937 global_rand;

        CliqueRunner(const InputGraph & g, const CliqueParams & p) :
            params(p),
            size(g.size()),
            adj(g.size(), BitSetType_{ unsigned(size), 0 }),
            order(size),
            invorder(size)
        {
            if (params.restarts_schedule->might_restart())
                watches.table.data.resize(g.size());

            // populate our order with every vertex initially
            iota(order.begin(), order.end(), 0);

            // pre-calculate degrees
            vector<int> degrees;
            degrees.resize(size);
            for (auto e = g.begin_edges(), e_end = g.end_edges() ; e != e_end ; ++e)
                ++degrees[e->first.first];

            // sort on degree
            sort(order.begin(), order.end(),
                    [&] (int a, int b) { return true ^ (degrees[a] < degrees[b] || (degrees[a] == degrees[b] && a > b)); });

            for (unsigned i = 0 ; i < order.size() ; ++i)
                invorder[order[i]] = i;

            for (auto e = g.begin_edges(), e_end = g.end_edges() ; e != e_end ; ++e)
                adj[invorder[e->first.first]].set(invorder[e->first.second]);
        }

        auto colour_class_order(
                const BitSetType_ & p,
                ArrayType_ & p_order,
                ArrayType_ & p_bounds,
                int & p_end) -> void
        {
            BitSetType_ p_left = p;      // not coloured yet
            unsigned colour = 0;         // current colour
            p_end = 0;

            // while we've things left to colour
            while (p_left.any()) {
                // next colour
                ++colour;
                // things that can still be given this colour
                BitSetType_ q = p_left;

                // while we can still give something this colour
                while (q.any()) {
                    // first thing we can colour
                    int v = q.find_first();
                    p_left.reset(v);
                    q.reset(v);

                    // can't give anything adjacent to this the same colour
                    q &= ~adj[v];

                    // record in result
                    p_bounds[p_end] = colour;
                    p_order[p_end] = v;
                    ++p_end;
                }
            }
        }

        auto post_nogood(
                const vector<int> & c)
        {
            Nogood<int> nogood;
            nogood.literals.assign(c.begin(), c.end());
            watches.post_nogood(move(nogood));
        }

        auto expand(
                unsigned long long & nodes,
                unsigned long long & find_nodes,
                unsigned long long & prove_nodes,
                vector<int> & c,
                BitSetType_ & p) -> SearchResult
        {
            ++nodes;
            ++prove_nodes;

            // initial colouring
            ArrayType_ p_order;
            ArrayType_ p_bounds;
            if constexpr (is_same<ArrayType_, vector<int> >::value) {
                p_order.resize(size);
                p_bounds.resize(size);
            }

            int p_end = 0;
            colour_class_order(p, p_order, p_bounds, p_end);

            // for each v in p... (v comes later)
            for (int n = p_end - 1 ; n >= 0 ; --n) {
                // bound, timeout or early exit?
                if (params.timeout->should_abort())
                    return SearchResult::Aborted;

                if (c.size() + p_bounds[n] <= incumbent.value)
                    break;

                // if we've used k colours to colour k vertices, it's a clique
                if (p_bounds[n] == n + 1) {
                    auto c_save = c;
                    for ( ; n >= 0 ; --n)
                        c.push_back(p_order[n]);
                    incumbent.update(c, find_nodes, prove_nodes);
                    c = move(c_save);

                    if (params.decide && incumbent.value >= *params.decide)
                        return SearchResult::DecidedTrue;

                    break;
                }

                auto v = p_order[n];

                // consider taking v
                c.push_back(v);

                if (params.decide) {
                    incumbent.update(c, find_nodes, prove_nodes);
                    if (incumbent.value >= *params.decide)
                        return SearchResult::DecidedTrue;
                }

                // filter p to contain vertices adjacent to v
                BitSetType_ new_p = p;
                new_p &= adj[v];

                if (params.restarts_schedule->might_restart())
                    watches.propagate(v,
                            [&] (int literal) { return c.end() == find(c.begin(), c.end(), literal); },
                            [&] (int literal) { new_p.reset(literal); }
                            );

                if (new_p.any()) {
                    switch (expand(nodes, find_nodes, prove_nodes, c, new_p)) {
                        case SearchResult::Aborted:
                            return SearchResult::Aborted;

                        case SearchResult::DecidedTrue:
                            return SearchResult::DecidedTrue;

                        case SearchResult::Complete:
                            break;

                        case SearchResult::Restart:
                            // restore assignments before posting nogoods, it's easier
                            c.pop_back();

                            // post nogoods for everything we've done so far
                            for (int m = p_end - 1 ; m > n ; --m) {
                                c.push_back(p_order[m]);
                                post_nogood(c);
                                c.pop_back();
                            }

                            return SearchResult::Restart;
                    }
                }
                else
                    incumbent.update(c, find_nodes, prove_nodes);

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

        auto run() -> CliqueResult
        {
            CliqueResult result;

            if (params.decide)
                incumbent.value = *params.decide - 1;

            // do the search
            bool done = false;
            unsigned number_of_restarts = 0;

            BitSetType_ p{ unsigned(size), 0 };
            for (int i = 0 ; i < size ; ++i)
                p.set(i);

            while (! done) {
                ++number_of_restarts;

                // start watching new nogoods
                done = watches.apply_new_nogoods(
                        [&] (int literal) { p.reset(literal); }
                        );

                if (done)
                    break;

                watches.clear_new_nogoods();

                auto new_p = p;
                vector<int> c;
                switch (expand(result.nodes, result.find_nodes, result.prove_nodes, c, new_p)) {
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

            result.clique.clear();
            for (auto & v : incumbent.c)
                result.clique.insert(order[v]);

            return result;
        }
    };
}

auto solve_clique_problem(const InputGraph & graph, const CliqueParams & params) -> CliqueResult
{
    return select_graph_size<CliqueRunner, CliqueResult>(AllGraphSizes(), graph, params);
}

