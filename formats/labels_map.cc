/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "formats/labels_map.hh"

auto LabelsMap::remap_vertex_label(int v) -> int
{
    auto r = _vertex_mapping.try_emplace(v, _next_vertex_label);
    if (r.second)
        ++_next_vertex_label;
    return r.first->second;
}

auto LabelsMap::remap_edge_label(int v) -> int
{
    auto r = _edge_mapping.try_emplace(v, _next_edge_label);
    if (r.second)
        ++_next_edge_label;
    return r.first->second;
}

