/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "graph.hh"

#include <algorithm>

using std::count_if;

Graph::Graph(int size)
{
    if (0 != size)
        resize(size);
}

auto Graph::_position(int a, int b) const -> AdjacencyMatrix::size_type
{
    return (a * _size) + b;
}

auto Graph::resize(int size) -> void
{
    _size = size;
    _adjacency.resize(size * size);
}

auto Graph::add_edge(int a, int b) -> void
{
    _adjacency[_position(a, b)] = 1;
    _adjacency[_position(b, a)] = 1;
}

auto Graph::adjacent(int a, int b) const -> bool
{
    return _adjacency[_position(a, b)];
}

auto Graph::size() const -> int
{
    return _size;
}

auto Graph::degree(int a) const -> int
{
    return count_if(
            _adjacency.begin() + a * _size,
            _adjacency.begin() + (a + 1) * _size,
            [] (int x) { return !!x; });
}

