#include <gss/formats/graph_file_error.hh>
#include <gss/formats/input_graph.hh>
#include <gss/utils/vertex_name_map.hh>

#include <algorithm>
#include <iterator>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using std::back_inserter;
using std::count_if;
using std::distance;
using std::find;
using std::function;
using std::isgraph;
using std::logic_error;
using std::make_optional;
using std::make_pair;
using std::map;
using std::max;
using std::move;
using std::nullopt;
using std::numeric_limits;
using std::optional;
using std::pair;
using std::string;
using std::string_view;
using std::to_string;
using std::transform;
using std::vector;

namespace
{
    auto sanity_check_name(string_view name, const char * const explanation) -> void
    {
        if (0 != count_if(name.begin(), name.end(), [](unsigned char c) { return ! isgraph(c); })) {
            string safe_name;
            transform(name.begin(), name.end(), back_inserter(safe_name), [](unsigned char c) { return isgraph(c) ? c : '?'; });
            throw GraphFileError("Suspicious input detected: " + string(explanation) + " '" + string(safe_name) + "' contains non-printable characters");
        }
    }
}

struct InputGraph::Imp
{
    struct Edge
    {
        string label;
        optional<long long> cost;
    };

    int size = 0;
    InputGraphProperties properties;

    // Keyed on the endpoint pair, so adjacent() and degree() do not care whether the
    // graph is a multigraph. A simple graph never has more than one Edge in a list; a
    // multigraph has one per label, since its edges are unique by (from, to, label).
    map<pair<int, int>, vector<Edge>> edges;
    vector<string> vertex_labels;
    vector<optional<long long>> vertex_costs;
    BiMap vertex_names;
    bool loopy = false;

    // replace says what re-adding an existing edge of a simple graph does: the labelled
    // overloads replace it, and the unlabelled add_edge(a, b) has always left it alone.
    auto add_one(int a, int b, string_view label, optional<long long> cost, bool replace = true) -> void
    {
        auto & list = edges[{a, b}];
        if (properties.multigraph) {
            for (auto & e : list)
                if (e.label == label) {
                    e.cost = cost;
                    return;
                }
            list.push_back(Edge{string{label}, cost});
        }
        else if (list.empty() || replace) {
            // A simple graph keeps one edge per pair.
            list.clear();
            list.push_back(Edge{string{label}, cost});
        }

        if (a == b)
            loopy = true;
    }

    auto check_cost_declaration(bool given) const -> void
    {
        if (given && ! properties.has_edge_costs)
            throw logic_error{"an edge was given a cost, but the graph was not declared to have edge costs"};
        if ((! given) && properties.has_edge_costs)
            throw logic_error{"an edge was added without a cost, but the graph was declared to have edge costs: an absent cost is not the same as 0"};
    }

    auto check_directed_declaration() const -> void
    {
        // Directedness is declared, not acquired. Adding a one-way edge to a graph that
        // says it is undirected would leave an asymmetric edge set behind a directed()
        // of false, which is the kind of quiet disagreement that shows up far from here.
        if (! properties.directed)
            throw logic_error{"add_directed_edge() on a graph that was not declared directed: pass directed = true to the InputGraph constructor, or use add_edge()"};
    }
};

InputGraph::InputGraph(int size, bool v, bool e, bool d) :
    InputGraph(size, InputGraphProperties{.has_vertex_labels = v, .has_edge_labels = e, .directed = d})
{
}

InputGraph::InputGraph(int size, const InputGraphProperties & properties) :
    _imp(std::make_unique<Imp>())
{
    _imp->properties = properties;

    if (0 != size)
        resize(size);
}

InputGraph::~InputGraph() = default;

InputGraph::InputGraph(InputGraph && other) = default;

auto InputGraph::resize(int size) -> void
{
    _imp->size = size;
    _imp->vertex_labels.resize(size);
    _imp->vertex_costs.resize(size);
}

auto InputGraph::add_edge(int a, int b) -> void
{
    _imp->check_cost_declaration(false);
    _imp->add_one(a, b, "", nullopt, false);
    _imp->add_one(b, a, "", nullopt, false);
}

auto InputGraph::add_edge(int a, int b, string_view label) -> void
{
    sanity_check_name(label, "edge label");
    _imp->check_cost_declaration(false);
    _imp->add_one(a, b, label, nullopt);
    _imp->add_one(b, a, label, nullopt);
}

auto InputGraph::add_directed_edge(int a, int b, string_view label) -> void
{
    sanity_check_name(label, "edge label");
    _imp->check_directed_declaration();
    _imp->check_cost_declaration(false);
    _imp->add_one(a, b, label, nullopt);
}

auto InputGraph::add_edge(int a, int b, string_view label, long long cost) -> void
{
    sanity_check_name(label, "edge label");
    _imp->check_cost_declaration(true);
    _imp->add_one(a, b, label, cost);
    _imp->add_one(b, a, label, cost);
}

auto InputGraph::add_directed_edge(int a, int b, string_view label, long long cost) -> void
{
    sanity_check_name(label, "edge label");
    _imp->check_directed_declaration();
    _imp->check_cost_declaration(true);
    _imp->add_one(a, b, label, cost);
}

auto InputGraph::adjacent(int a, int b) const -> bool
{
    return _imp->edges.count({a, b});
}

auto InputGraph::size() const -> int
{
    return _imp->size;
}

auto InputGraph::number_of_directed_edges() const -> int
{
    int result = 0;
    for (auto & [_, list] : _imp->edges)
        result += list.size();
    return result;
}

auto InputGraph::loopy() const -> bool
{
    return _imp->loopy;
}

auto InputGraph::degree(int a) const -> int
{
    auto lower = _imp->edges.lower_bound({a, 0});
    auto upper = _imp->edges.upper_bound({a, numeric_limits<int>::max()});
    return distance(lower, upper);
}

auto InputGraph::set_vertex_label(int v, string_view l) -> void
{
    sanity_check_name(l, "vertex label");
    _imp->vertex_labels[v] = l;
}

auto InputGraph::vertex_label(int v) const -> string_view
{
    return _imp->vertex_labels[v];
}

auto InputGraph::set_vertex_name(int v, string_view l) -> void
{
    sanity_check_name(l, "vertex name");
    _imp->vertex_names.erase(v);
    _imp->vertex_names.insert(v, string(l));
}

auto InputGraph::vertex_name(int v) const -> string
{
    auto it = _imp->vertex_names.find_left(v);
    if (it == _imp->vertex_names.id_to_name.end())
        return to_string(v);
    else
        return it->second;
}

auto InputGraph::vertex_has_name(int v) const -> bool
{
    return _imp->vertex_names.find_left(v) != _imp->vertex_names.id_to_name.end();
}

auto InputGraph::vertex_from_name(string_view n) const -> optional<int>
{
    auto it = _imp->vertex_names.find_right(string(n));
    if (it == _imp->vertex_names.name_to_id.end())
        return nullopt;
    else
        return make_optional(it->second);
}

auto InputGraph::edge_label(int a, int b) const -> string_view
{
    if (_imp->properties.multigraph)
        throw logic_error{"edge_label() on a multigraph, where a pair of vertices may have several edges: use for_each_edge()"};
    return _imp->edges.find({a, b})->second.front().label;
}

auto InputGraph::has_vertex_labels() const -> bool
{
    return _imp->properties.has_vertex_labels;
}

auto InputGraph::has_edge_labels() const -> bool
{
    return _imp->properties.has_edge_labels;
}

auto InputGraph::directed() const -> bool
{
    return _imp->properties.directed;
}

auto InputGraph::multigraph() const -> bool
{
    return _imp->properties.multigraph;
}

auto InputGraph::has_vertex_costs() const -> bool
{
    return _imp->properties.has_vertex_costs;
}

auto InputGraph::has_edge_costs() const -> bool
{
    return _imp->properties.has_edge_costs;
}

auto InputGraph::set_vertex_cost(int v, long long cost) -> void
{
    if (! _imp->properties.has_vertex_costs)
        throw logic_error{"set_vertex_cost() on a graph that was not declared to have vertex costs"};
    _imp->vertex_costs[v] = cost;
}

auto InputGraph::vertex_cost(int v) const -> long long
{
    if (! _imp->properties.has_vertex_costs)
        throw logic_error{"vertex_cost() on a graph that was not declared to have vertex costs"};
    if (! _imp->vertex_costs[v])
        throw logic_error{"vertex " + to_string(v) + " was never given a cost, but the graph was declared to have vertex costs: an absent cost is not the same as 0"};
    return *_imp->vertex_costs[v];
}

auto InputGraph::for_each_edge(const function<auto(int, int, std::string_view)->void> & c) const -> void
{
    for (auto & [e, list] : _imp->edges)
        for (auto & edge : list)
            c(e.first, e.second, edge.label);
}

auto InputGraph::for_each_edge_and_cost(const function<auto(int, int, std::string_view, optional<long long>)->void> & c) const -> void
{
    for (auto & [e, list] : _imp->edges)
        for (auto & edge : list)
            c(e.first, e.second, edge.label, edge.cost);
}
