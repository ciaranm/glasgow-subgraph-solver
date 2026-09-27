#include <gss/innards/reification.hh>

#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

using namespace gss::innards;

using std::optional;
using std::string;
using std::string_view;
using std::tuple;
using std::vector;

namespace
{
    struct Edge
    {
        int from, to;
        string label;
        long long cost;
    };

    // Every edge of g that becomes an edge-vertex. In an undirected reification an
    // undirected edge is one edge-vertex, so only one of its two directions is kept; in
    // a directed one it is two arcs and so two edge-vertices. A loop is one either way.
    auto edges_to_reify(const InputGraph & g, bool directed) -> vector<Edge>
    {
        vector<Edge> result;
        g.for_each_edge_and_cost([&](int f, int t, string_view l, optional<long long> c) {
            if ((! directed) && t < f)
                return;
            result.push_back(Edge{f, t, string{l}, c.value_or(0)});
        });
        return result;
    }

    auto reify_one(const InputGraph & g, bool directed, bool use_vertex_labels, bool use_edge_labels,
        vector<EdgeVertexEndpoints> & endpoints, vector<long long> * costs) -> InputGraph
    {
        auto edges = edges_to_reify(g, directed);
        int n = g.size();

        InputGraph result{int(n + edges.size()), InputGraphProperties{.has_vertex_labels = true, .directed = directed}};

        for (int v = 0; v < n; ++v) {
            result.set_vertex_label(v, "v" + (use_vertex_labels ? string{g.vertex_label(v)} : string{}));
            if (g.vertex_has_name(v))
                result.set_vertex_name(v, g.vertex_name(v));
            if (costs)
                costs->push_back(g.has_vertex_costs() ? g.vertex_cost(v) : 0);
        }

        for (size_t i = 0; i < edges.size(); ++i) {
            auto & [f, t, l, c] = edges[i];
            int e = n + int(i);
            result.set_vertex_label(e, (f == t ? "l" : "e") + (use_edge_labels ? l : string{}));
            if (directed) {
                result.add_directed_edge(f, e, "");
                result.add_directed_edge(e, t, "");
            }
            else {
                result.add_edge(f, e);
                result.add_edge(e, t);
            }
            endpoints.push_back(EdgeVertexEndpoints{f, t});
            if (costs)
                costs->push_back(c);
        }

        return result;
    }
}

auto gss::innards::needs_reification(const InputGraph & pattern, const InputGraph & target, bool minimising_cost) -> bool
{
    return pattern.multigraph() || target.multigraph() || (minimising_cost && target.has_edge_costs());
}

auto gss::innards::reify(const InputGraph & pattern, const InputGraph & target) -> Reification
{
    bool directed = pattern.directed() || target.directed();

    // The model ignores the target's labels of a kind the pattern does not have, and so
    // does this, so that reifying changes nothing about which mappings are allowed.
    bool use_vertex_labels = pattern.has_vertex_labels();
    bool use_edge_labels = pattern.has_edge_labels();

    vector<EdgeVertexEndpoints> pattern_edge_vertices, target_edge_vertices;
    vector<long long> target_costs;

    auto reified_pattern = reify_one(pattern, directed, use_vertex_labels, use_edge_labels, pattern_edge_vertices, nullptr);
    auto reified_target = reify_one(target, directed, use_vertex_labels, use_edge_labels, target_edge_vertices, &target_costs);

    return Reification{
        std::move(reified_pattern),
        std::move(reified_target),
        pattern.size(),
        target.size(),
        std::move(pattern_edge_vertices),
        std::move(target_edge_vertices),
        directed,
        std::move(target_costs)};
}
