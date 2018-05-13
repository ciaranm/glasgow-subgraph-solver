/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_FORMATS_LABELS_MAP_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_FORMATS_LABELS_MAP_HH 1

#include <map>

class LabelsMap
{
    private:
        std::map<int, int> _vertex_mapping, _edge_mapping;
        int _next_vertex_label = 0, _next_edge_label = 0;

    public:
        auto remap_vertex_label(int l) -> int;
        auto remap_edge_label(int l) -> int;
};

#endif
