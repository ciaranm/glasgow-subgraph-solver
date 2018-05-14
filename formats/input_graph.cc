/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "formats/input_graph.hh"

#include <algorithm>
#include <limits>
#include <string_view>

using std::distance;
using std::numeric_limits;
using std::make_pair;
using std::max;
using std::string_view;

InputGraph::InputGraph(int size, bool v, bool e) :
    _has_vertex_labels(v),
    _has_edge_labels(e)
{
    if (0 != size)
        resize(size);
}

auto InputGraph::resize(int size) -> void
{
    _size = size;
    _vertex_labels.resize(size);
}

auto InputGraph::add_edge(int a, int b) -> void
{
    _edges.emplace(make_pair(a, b), "");
    _edges.emplace(make_pair(b, a), "");
}

auto InputGraph::add_directed_edge(int a, int b, string_view label) -> void
{
    _edges.emplace(make_pair(a, b), label).first->second = label;
    _edges.emplace(make_pair(b, a), "unlabelled");
}

auto InputGraph::adjacent(int a, int b) const -> bool
{
    return _edges.count({ a, b });
}

auto InputGraph::size() const -> int
{
    return _size;
}

auto InputGraph::degree(int a) const -> int
{
    auto lower = _edges.lower_bound({ a, 0 });
    auto upper = _edges.upper_bound({ a, numeric_limits<int>::max() });
    return distance(lower, upper);
}

auto InputGraph::set_vertex_label(int v, string_view l) -> void
{
    _vertex_labels[v] = l;
}

auto InputGraph::vertex_label(int v) const -> string_view
{
    return _vertex_labels[v];
}

auto InputGraph::edge_label(int a, int b) const -> string_view
{
    return _edges.find({a, b})->second;
}

auto InputGraph::begin_edges() const -> InputGraph::EdgesIterator
{
    return _edges.begin();
}

auto InputGraph::end_edges() const -> InputGraph::EdgesIterator
{
    return _edges.end();
}

auto InputGraph::has_vertex_labels() const -> bool
{
    return _has_vertex_labels;
}

auto InputGraph::has_edge_labels() const -> bool
{
    return _has_edge_labels;
}

