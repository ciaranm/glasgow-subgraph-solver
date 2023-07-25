#include <gss/formats/graph_file_error.hh>
#include <gss/formats/input_graph.hh>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <iterator>
#include <limits>
#include <map>
#include <string>
#include <type_traits>
#include <vector>

#include <boost/bimap.hpp>
#include <boost/bimap/unordered_set_of.hpp>
#include <boost/container/allocator.hpp>

using std::back_inserter;
using std::count_if;
using std::distance;
using std::find;
using std::function;
using std::isgraph;
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

using Names = boost::bimap<
    boost::bimaps::unordered_set_of<int>,
    boost::bimaps::unordered_set_of<string>,
    boost::container::allocator<pair<int, string>>>;

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
    int size = 0;
    bool has_vertex_labels, has_edge_labels;
    map<pair<int, int>, string> edges;
    vector<string> vertex_labels;
    Names vertex_names;
    bool loopy = false, directed = false;
};

InputGraph::InputGraph(int size, bool v, bool e) :
    _imp(new Imp{})
{
    _imp->has_vertex_labels = v;
    _imp->has_edge_labels = e;
    _imp->loopy = false;
    _imp->directed = false;

    if (0 != size)
        resize(size);
}

InputGraph::~InputGraph() = default;

InputGraph::InputGraph(InputGraph && other) = default;

auto InputGraph::resize(int size) -> void
{
    _imp->size = size;
    _imp->vertex_labels.resize(size);
}

auto InputGraph::add_edge(int a, int b) -> void
{
    _imp->edges.emplace(make_pair(a, b), "");
    _imp->edges.emplace(make_pair(b, a), "");
    if (a == b)
        _imp->loopy = true;
}

auto InputGraph::add_directed_edge(int a, int b, string_view label) -> void
{
    sanity_check_name(label, "edge label");

    _imp->directed = true;

    _imp->edges.emplace(make_pair(a, b), label).first->second = label;
    if (a == b)
        _imp->loopy = true;
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
    return _imp->edges.size();
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
    _imp->vertex_names.left.erase(v);
    _imp->vertex_names.insert(Names::value_type{v, string{l}});
}

auto InputGraph::vertex_name(int v) const -> string
{
    auto it = _imp->vertex_names.left.find(v);
    if (it == _imp->vertex_names.left.end())
        return to_string(v);
    else
        return it->second;
}

auto InputGraph::vertex_from_name(string_view n) const -> optional<int>
{
    auto it = _imp->vertex_names.right.find(string{n});
    if (it == _imp->vertex_names.right.end())
        return nullopt;
    else
        return make_optional(it->second);
}

auto InputGraph::edge_label(int a, int b) const -> string_view
{
    return _imp->edges.find({a, b})->second;
}

auto InputGraph::has_vertex_labels() const -> bool
{
    return _imp->has_vertex_labels;
}

auto InputGraph::has_edge_labels() const -> bool
{
    return _imp->has_edge_labels;
}

auto InputGraph::directed() const -> bool
{
    return _imp->directed;
}

auto InputGraph::for_each_edge(const function<auto(int, int, std::string_view)->void> & c) const -> void
{
    for (auto & [e, l] : _imp->edges)
        c(e.first, e.second, l);
}
