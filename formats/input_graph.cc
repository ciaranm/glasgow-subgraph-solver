/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "formats/input_graph.hh"

#include <algorithm>
#include <limits>

using std::distance;
using std::numeric_limits;

InputGraph::InputGraph(int size)
{
    if (0 != size)
        resize(size);
}

auto InputGraph::resize(int size) -> void
{
    _size = size;
}

auto InputGraph::add_edge(int a, int b) -> void
{
    _edges.emplace(a, b);
    _edges.emplace(b, a);
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

