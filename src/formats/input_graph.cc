/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "formats/input_graph.hh"

#include <algorithm>
#include <limits>
#include <string_view>
#include <vector>

using std::distance;
using std::find;
using std::numeric_limits;
using std::make_optional;
using std::make_pair;
using std::max;
using std::nullopt;
using std::optional;
using std::pair;
using std::string;
using std::string_view;
using std::to_string;
using std::vector;

InputGraph::InputGraph(int size, bool v, bool e, bool f) :
    _has_vertex_labels(v),
    _has_edge_labels(e),
    _loopy(false),
    _directed(f)
{
    if (0 != size)
        resize(size);
}

auto InputGraph::resize(int size) -> void
{
    _size = size;
    _vertex_labels.resize(size);
    _vertex_names.resize(size);
    _vertex_pattern_constraints.resize(size);
    _vertex_directed_degrees.resize(size);
}

auto InputGraph::add_edge(int a, int b) -> void
{
    _edges.emplace(make_pair(a, b), "");
    _edges.emplace(make_pair(b, a), "");
    if (a == b)
        _loopy = true;
}

auto InputGraph::add_directed_edge(int a, int b, const string & label) -> void
{
    _directed = true;
    if (a == b)
        _loopy = true;

    _edges.emplace(make_pair(a, b), label).first->second = label;

    if(vertex_label(a) != "LINK" && vertex_label(b) != "LINK"){ 
        _vertex_directed_degrees[b].first++;
        _vertex_directed_degrees[a].second++;
    }
}

auto InputGraph::adjacent(int a, int b) const -> bool
{
    return _edges.count({ a, b });
}

auto InputGraph::add_link_adjacency(int a, int b, int c) -> void
{
    if(c < 1) return;

    if(_link_adjacencies[_no_link_nodes+b-_size].find(_vertex_labels[a]) == _link_adjacencies[_no_link_nodes+b-_size].end())
        _link_adjacencies[_no_link_nodes+b-_size][_vertex_labels[a]] = vector<int> {};

    _link_adjacencies[_no_link_nodes+b-_size].find(_vertex_labels[a])->second.push_back(c);
}

auto InputGraph::get_link_adjacency_list(int a) const -> std::map<std::string, std::vector<int>>
{
    return _link_adjacencies[_no_link_nodes+a-_size];
}

auto InputGraph::size() const -> int
{
    return _size;
}

auto InputGraph::number_of_directed_edges() const -> int
{
    return _edges.size();
}

auto InputGraph::add_pattern_site_edge(int a, int b) -> void
{
    _pattern_site_edges.push_back(make_pair(a, b));
}

auto InputGraph::get_pattern_site_edge(int s) const -> pair<int, int>
{
    return _pattern_site_edges[s];
}

auto InputGraph::no_pattern_site_edges() const -> int
{
    return _pattern_site_edges.size();
}


auto InputGraph::add_pattern_root_edge(int a, int b) -> void
{
    _pattern_root_edges.push_back(make_pair(a, b));
}

auto InputGraph::get_pattern_root_edge(int r) const -> pair<int, int>
{
    return _pattern_root_edges[r];
}

auto InputGraph::no_pattern_root_edges() const -> int
{
    return _pattern_root_edges.size();
}

auto InputGraph::add_link_node() -> void
{
    resize(_size+1);
    set_vertex_label(_size-1, "LINK");
    //_link_adjacencies.resize(_link_adjacencies.size()+1);
    _no_link_nodes++;
}

auto InputGraph::get_no_link_nodes() const -> int
{
    return _no_link_nodes;
}

auto InputGraph::loopy() const -> bool
{
    return _loopy;
}

auto InputGraph::degree(int a) const -> int
{
    auto lower = _edges.lower_bound({ a, 0 });
    auto upper = _edges.upper_bound({ a, numeric_limits<int>::max() });
    return distance(lower, upper);
}

auto InputGraph::in_degree(int a) const -> int
{
    return _vertex_directed_degrees[a].first;
}

auto InputGraph::out_degree(int a) const -> int
{
    return _vertex_directed_degrees[a].second;
}

auto InputGraph::set_vertex_label(int v, string_view l) -> void
{
    _vertex_labels[v] = l;
}

auto InputGraph::vertex_label(int v) const -> string_view
{
    return _vertex_labels[v];
}

auto InputGraph::set_vertex_name(int v, string_view l) -> void
{
    _vertex_names[v] = l;
}

auto InputGraph::set_child_of_root(int v) -> void
{
    _vertex_pattern_constraints.at(v).first = true;
}

auto InputGraph::set_parent_of_site(int v) -> void
{
    _vertex_pattern_constraints.at(v).second = true;
}

auto InputGraph::get_big_constraint(int v) const -> pair<bool, bool>
{
    return _vertex_pattern_constraints.at(v);
}

auto InputGraph::vertex_name(int v) const -> string
{
    if (_vertex_names[v].empty())
        return to_string(v);
    else
        return _vertex_names[v];
}

auto InputGraph::vertex_from_name(string_view n) const -> optional<int>
{
    auto i = find(_vertex_names.begin(), _vertex_names.end(), n);
    if (i == _vertex_names.end())
        return nullopt;
    else
        return make_optional(i - _vertex_names.begin());
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

auto InputGraph::directed() const -> bool
{
    return _directed;
}

