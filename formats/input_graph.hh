/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_SOLVER_FORMATS_INPUT_GRAPH_HH
#define GLASGOW_SUBGRAPH_SOLVER_SOLVER_FORMATS_INPUT_GRAPH_HH 1

#include <cstdint>
#include <map>
#include <string>
#include <type_traits>
#include <vector>

/**
 * A graph, in a convenient format for reading in from files. We don't do any
 * performance critical operations on this: the algorithms re-encode as
 * necessary.
 *
 * Indices start at 0.
 */
class InputGraph
{
    private:
        int _size = 0;
        bool _has_vertex_labels, _has_edge_labels;
        std::map<std::pair<int, int>, std::string> _edges;
        std::vector<std::string> _vertex_labels;

    public:
        /**
         * \param initial_size can be 0, if resize() is called afterwards.
         */
        InputGraph(int initial_size, bool has_vertex_labels, bool has_edge_labels);

        InputGraph(const InputGraph &) = default;

        /**
         * Number of vertices.
         */
        auto size() const -> int;

        /**
         * Change our size. Must be called before adding an edge.
         */
        auto resize(int size) -> void;

        /**
         * Add an edge from a to b (and from b to a).
         */
        auto add_edge(int a, int b) -> void;

        /**
         * Add a directed edge from a to b, with a label.
         */
        auto add_directed_edge(int a, int b, std::string_view label) -> void;

        /**
         * Are vertices a and b adjacent?
         */
        auto adjacent(int a, int b) const -> bool;

        /**
         * What is the degree of a given vertex?
         */
        auto degree(int a) const -> int;

        /**
         * Set a vertex label.
         */
        auto set_vertex_label(int v, std::string_view label) -> void;

        /**
         * What is the label associated with a given vertex?
         */
        auto vertex_label(int v) const -> std::string_view;

        auto has_vertex_labels() const -> bool;

        /**
         * What is the label associated with a given edge?
         */
        auto edge_label(int a, int b) const -> std::string_view;

        auto has_edge_labels() const -> bool;

        using EdgesIterator = std::map<std::pair<int, int>, std::string>::const_iterator;

        auto begin_edges() const -> EdgesIterator;

        auto end_edges() const -> EdgesIterator;
};

#endif
