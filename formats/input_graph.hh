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
        int _last_vertex_label = 0, _last_edge_label = 0;
        std::map<std::pair<int, int>, int> _edges;
        std::vector<int> _vertex_labels;

    public:
        /**
         * \param initial_size can be 0, if resize() is called afterwards.
         */
        InputGraph(int initial_size);

        InputGraph(const InputGraph &) = default;

        explicit InputGraph() = default;

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
        auto add_edge(int a, int b, int label = 0) -> void;

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
        auto set_vertex_label(int v, int l) -> void;

        /**
         * What is the label associated with a given vertex?
         *
         * In the range 0..vertex_labels_count() - 1, always 0 for unlabelled graphs.
         */
        auto vertex_label(int v) const -> int;

        /**
         * How many vertex labels do we have? For unlabelled graphs, every vertex
         * has label 0, and so we have 1.
         */
        auto vertex_labels_count() const -> int;

        /**
         * What is the label associated with a given edge?
         *
         * In the range 0..edge_labels_count() - 1, always 0 for unlabelled graphs.
         */
        auto edge_label(int a, int b) const -> int;

        /**
         * How many edge labels do we have? For unlabelled graphs, every edge
         * has label 0, and so we have 1.
         */
        auto edge_labels_count() const -> int;
};

#endif
