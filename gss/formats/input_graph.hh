#ifndef GLASGOW_SUBGRAPH_SOLVER_SOLVER_FORMATS_INPUT_GRAPH_HH
#define GLASGOW_SUBGRAPH_SOLVER_SOLVER_FORMATS_INPUT_GRAPH_HH 1

#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <string_view>

/**
 * The graph-level properties of an InputGraph. Every one of these is declared, never
 * inferred from which vertices, edges, labels or costs happen to get added: see
 * dev_docs/file-formats.md for why. Each defaults to the reading under which a caller
 * that should have declared it gets an error rather than a graph that quietly means
 * something else.
 */
struct InputGraphProperties
{
    bool has_vertex_labels = false;
    bool has_edge_labels = false;
    bool directed = false;

    /**
     * May a pair of vertices be joined by several edges, with different labels? Edges
     * are then unique by (from, to, label) rather than by (from, to), and edge_label()
     * refuses, since there is no single answer to give. Only the homomorphism solver
     * accepts a multigraph, and it rewrites one before search (see
     * gss/innards/reification.hh).
     */
    bool multigraph = false;

    /**
     * Does every vertex carry an integer cost? If so every vertex must be given one,
     * with set_vertex_cost(): an absent cost is not the same as 0.
     */
    bool has_vertex_costs = false;

    /**
     * Does every edge carry an integer cost? If so every edge must be added with one of
     * the overloads that takes a cost, and the others refuse.
     */
    bool has_edge_costs = false;
};

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
    struct Imp;
    std::unique_ptr<Imp> _imp;

public:
    /**
     * \param initial_size can be 0, if resize() is called afterwards.
     * \param directed declares whether this is a directed graph. It is never
     *     inferred from which edges get added: add_directed_edge() requires it to
     *     have been declared, and an undirected graph stays undirected however its
     *     edges are labelled. Defaults to false, the reading under which a caller
     *     that should have declared it gets an error rather than a graph that
     *     quietly means something else.
     */
    InputGraph(int initial_size, bool has_vertex_labels, bool has_edge_labels, bool directed = false);

    /**
     * \param initial_size can be 0, if resize() is called afterwards.
     */
    InputGraph(int initial_size, const InputGraphProperties & properties);

    InputGraph(const InputGraph &) = delete;

    InputGraph(InputGraph &&);

    ~InputGraph();

    /**
     * Number of vertices.
     */
    auto size() const -> int;

    /**
     * Number of (directed, even if the graph is undirected) edges.
     */
    auto number_of_directed_edges() const -> int;

    /**
     * Do we have any loops?
     */
    auto loopy() const -> bool;

    /**
     * Change our size. Must be called before adding an edge.
     */
    auto resize(int size) -> void;

    /**
     * Add an edge from a to b (and from b to a).
     */
    auto add_edge(int a, int b) -> void;

    /**
     * Add an edge from a to b (and from b to a), with a label.
     *
     * Unlike add_directed_edge(), this does not make the graph directed: an
     * undirected graph whose edges happen to be labelled is still undirected.
     */
    auto add_edge(int a, int b, std::string_view label) -> void;

    /**
     * Add a directed edge from a to b, with a label.
     *
     * \throw std::logic_error if the graph was not declared directed.
     */
    auto add_directed_edge(int a, int b, std::string_view label) -> void;

    /**
     * Add an edge from a to b (and from b to a), with a label and a cost.
     *
     * \throw std::logic_error if the graph was not declared to have edge costs.
     */
    auto add_edge(int a, int b, std::string_view label, long long cost) -> void;

    /**
     * Add a directed edge from a to b, with a label and a cost.
     *
     * \throw std::logic_error if the graph was not declared directed, or not declared
     *     to have edge costs.
     */
    auto add_directed_edge(int a, int b, std::string_view label, long long cost) -> void;

    /**
     * Are vertices a and b adjacent?
     */
    auto adjacent(int a, int b) const -> bool;

    /**
     * What is the degree of a given vertex? In a multigraph this counts neighbours, not
     * edges: parallel edges to one neighbour count once.
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
     * Set a vertex name (for output purposes).
     */
    auto set_vertex_name(int v, std::string_view label) -> void;

    /**
     * What is the name associated with a given vertex (for output purposes)?
     */
    auto vertex_name(int v) const -> std::string;

    /**
     * Has a name actually been set for this vertex? vertex_name() falls back to the
     * index when one hasn't, which a writer needs to tell apart from a vertex
     * genuinely named for its index.
     */
    auto vertex_has_name(int v) const -> bool;

    /**
     * Find a given vertex by name.
     */
    auto vertex_from_name(std::string_view n) const -> std::optional<int>;

    /**
     * What is the label associated with a given edge?
     *
     * \throw std::logic_error on a multigraph, where a pair of vertices may have several.
     */
    auto edge_label(int a, int b) const -> std::string_view;

    auto has_edge_labels() const -> bool;

    auto directed() const -> bool;

    auto multigraph() const -> bool;

    auto has_vertex_costs() const -> bool;

    auto has_edge_costs() const -> bool;

    /**
     * Set a vertex's cost.
     *
     * \throw std::logic_error if the graph was not declared to have vertex costs.
     */
    auto set_vertex_cost(int v, long long cost) -> void;

    /**
     * What is a vertex's cost?
     *
     * \throw std::logic_error if the graph was not declared to have vertex costs, or
     *     if this vertex was never given one.
     */
    auto vertex_cost(int v) const -> long long;

    /**
     * Every edge, including each of several parallel edges in a multigraph. An
     * undirected edge is visited once in each direction, and a loop once.
     */
    auto for_each_edge(const std::function<auto(int, int, std::string_view)->void> &) const -> void;

    /**
     * As for_each_edge(), with each edge's cost, or nullopt if the graph was not
     * declared to have edge costs.
     */
    auto for_each_edge_and_cost(const std::function<auto(int, int, std::string_view, std::optional<long long>)->void> &) const -> void;
};

#endif
