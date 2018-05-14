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

InputGraph::InputGraph(int size)
{
    if (0 != size)
        resize(size);
}

auto InputGraph::resize(int size) -> void
{
    _size = size;
    _vertex_labels.resize(size);
}

auto InputGraph::add_edge(int a, int b, string_view label) -> void
{
    _edges.emplace(make_pair(a, b), label);
    _edges.emplace(make_pair(b, a), label);
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

