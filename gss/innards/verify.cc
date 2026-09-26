#include <gss/innards/verify.hh>

#include <map>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <tuple>

using namespace gss;
using namespace gss::innards;

using std::map;
using std::optional;
using std::set;
using std::string;
using std::string_view;
using std::tuple;

namespace
{
    // Every target edge, with its label and cost, keyed so that a pattern edge can look
    // up the one it must land on. Works for multigraphs, where edge_label() does not.
    auto target_edges(const InputGraph & target) -> map<tuple<int, int, string>, optional<long long>>
    {
        map<tuple<int, int, string>, optional<long long>> result;
        target.for_each_edge_and_cost([&](int f, int t, string_view l, optional<long long> c) {
            result.emplace(tuple{f, t, string{l}}, c);
        });
        return result;
    }
}

BuggySolution::BuggySolution(const string & message) noexcept :
    _what(message)
{
}

auto BuggySolution::what() const noexcept -> const char *
{
    return _what.c_str();
}

auto gss::innards::verify_homomorphism(
    const InputGraph & pattern,
    const InputGraph & target,
    bool injective,
    bool locally_injective,
    bool induced,
    const map<int, int> & mapping) -> void
{
    // nothing to verify, if unsat
    if (mapping.empty())
        return;

    // totality and range
    for (int i = 0; i < pattern.size(); ++i) {
        if (! mapping.count(i))
            throw BuggySolution{"No mapping for vertex " + pattern.vertex_name(i)};
        else if (mapping.find(i)->second < 0 || mapping.find(i)->second >= target.size())
            throw BuggySolution{"Mapping " + pattern.vertex_name(i) + " -> " +
                target.vertex_name(mapping.find(i)->second) + " out of range"};
    }

    // no extra stuff
    for (auto & i : mapping)
        if (i.first < 0 || i.first >= pattern.size())
            throw BuggySolution{"Vertex " + pattern.vertex_name(i.first) + " out of range"};

    // labels
    if (pattern.has_vertex_labels())
        for (int i = 0; i < pattern.size(); ++i)
            if (pattern.vertex_label(i) != target.vertex_label(mapping.find(i)->second))
                throw BuggySolution{"Mismatched vertex label for assignment " + pattern.vertex_name(i) + " -> " +
                    target.vertex_name(mapping.find(i)->second)};

    // injectivity
    if (injective) {
        map<int, int> seen;
        for (auto & [i, j] : mapping) {
            if (! seen.emplace(j, i).second)
                throw BuggySolution{"Non-injective mapping: " + pattern.vertex_name(i) + " -> " +
                    target.vertex_name(mapping.find(i)->second) + " and " +
                    pattern.vertex_name(seen.find(j)->second) + " -> " + target.vertex_name(j)};
        }
    }

    // local injectivity
    if (locally_injective) {
        for (int v = 0; v < pattern.size(); ++v) {
            map<int, int> seen;
            for (auto & [i, j] : mapping) {
                if (pattern.adjacent(v, i) && ! seen.emplace(j, i).second)
                    throw BuggySolution{"Non locally-injective mapping: on neighbourhood of " + pattern.vertex_name(v) + ", " + pattern.vertex_name(i) + " -> " +
                        target.vertex_name(mapping.find(i)->second) + " and " +
                        pattern.vertex_name(seen.find(j)->second) + " -> " + target.vertex_name(j)};
            }
        }
    }

    // loops
    for (int i = 0; i < pattern.size(); ++i) {
        if (pattern.adjacent(i, i) && ! target.adjacent(mapping.find(i)->second, mapping.find(i)->second))
            throw BuggySolution{"Vertex " + pattern.vertex_name(i) + " has a loop but mapped vertex " +
                target.vertex_name(mapping.find(i)->second) + " does not"};
        else if (induced && target.adjacent(mapping.find(i)->second, mapping.find(i)->second) &&
            ! pattern.adjacent(i, i))
            throw BuggySolution{"Vertex " + pattern.vertex_name(i) + " has no loop but mapped vertex " +
                target.vertex_name(mapping.find(i)->second) + " does"};
    }

    // adjacency, non-adjacency, and edge labels
    for (auto & [i, t] : mapping) {
        for (auto & [j, u] : mapping) {
            if (pattern.adjacent(i, j) && ! target.adjacent(t, u))
                throw BuggySolution{"Edge " + pattern.vertex_name(i) + " -- " + pattern.vertex_name(j) +
                    " mapped to non-edge " + target.vertex_name(t) + " -/- " + target.vertex_name(u)};
            else if (induced && ! pattern.adjacent(i, j) && target.adjacent(t, u))
                throw BuggySolution{"Non-edge " + pattern.vertex_name(i) + " -/- " + pattern.vertex_name(j) +
                    " mapped to edge " + target.vertex_name(t) + " -- " + target.vertex_name(u)};
            // Edge labels were not checked here at all, despite the heading. The pairs
            // with i == j are in this cross product, so this covers a loop's label too --
            // which is the one the solver itself can miss, loops being stripped out of its
            // adjacency rows (issue #92). A multigraph's labels are checked below instead,
            // since a pair may have several.
            else if (pattern.has_edge_labels() && ! pattern.multigraph() && ! target.multigraph() && pattern.adjacent(i, j) &&
                pattern.edge_label(i, j) != target.edge_label(t, u))
                throw BuggySolution{"Edge " + pattern.vertex_name(i) + " -- " + pattern.vertex_name(j) +
                    " labelled '" + string{pattern.edge_label(i, j)} + "' mapped to edge " +
                    target.vertex_name(t) + " -- " + target.vertex_name(u) + " labelled '" +
                    string{target.edge_label(t, u)} + "'"};
        }
    }

    // In a multigraph, every labelled pattern edge needs a target edge with that label.
    if (pattern.has_edge_labels() && (pattern.multigraph() || target.multigraph())) {
        auto edges = target_edges(target);
        pattern.for_each_edge([&](int i, int j, string_view l) {
            auto t = mapping.find(i)->second, u = mapping.find(j)->second;
            if (! edges.contains(tuple{t, u, string{l}}))
                throw BuggySolution{"Edge " + pattern.vertex_name(i) + " -- " + pattern.vertex_name(j) +
                    " labelled '" + string{l} + "' mapped to " + target.vertex_name(t) + " -- " +
                    target.vertex_name(u) + ", which has no edge with that label"};
        });
    }
}

auto gss::innards::cost_of_mapping(
    const InputGraph & pattern,
    const InputGraph & target,
    const map<int, int> & mapping) -> long long
{
    long long result = 0;
    if (target.has_vertex_costs())
        for (auto & [_, t] : mapping)
            result += target.vertex_cost(t);

    if (target.has_edge_costs()) {
        bool directed = pattern.directed() || target.directed();
        auto edges = target_edges(target);
        bool use_labels = pattern.has_edge_labels();
        pattern.for_each_edge([&](int i, int j, string_view l) {
            if ((! directed) && j < i)
                return;
            auto t = mapping.find(i)->second, u = mapping.find(j)->second;
            // Without pattern edge labels, target labels are ignored, so any target edge
            // between the images will do, and the solver will have taken the cheapest.
            optional<long long> best;
            for (auto & [key, c] : edges)
                if (get<0>(key) == t && get<1>(key) == u && ((! use_labels) || get<2>(key) == l))
                    if ((! best) || *c < *best)
                        best = *c;
            if (! best)
                throw BuggySolution{"Edge " + pattern.vertex_name(i) + " -- " + pattern.vertex_name(j) +
                    " mapped to " + target.vertex_name(t) + " -- " + target.vertex_name(u) + ", which has no matching edge"};
            result += *best;
        });
    }

    return result;
}
