/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_COMMON_SUBGRAPH_HH
#define GLASGOW_SUBGRAPH_SOLVER_COMMON_SUBGRAPH_HH 1

#include "formats/input_graph.hh"
#include "timeout.hh"
#include "vertex_to_vertex_mapping.hh"
#include "proof-fwd.hh"
#include "loooong.hh"

#include <functional>
#include <list>
#include <memory>
#include <string>

struct CommonSubgraphParams
{
    /// Timeout handler
    std::shared_ptr<Timeout> timeout;

    /// The start time of the algorithm.
    std::chrono::time_point<std::chrono::steady_clock> start_time;

    /// Decide instead of maximise?
    std::optional<unsigned> decide;

    /// Enumerate? (Decide only)
    bool count_solutions = false;

    /// Print solutions, for enumerating
    std::function<auto (const VertexToVertexMapping &) -> void> enumerate_callback;

    /// Optional proof handler
    std::shared_ptr<Proof> proof;

    /// Connected?
    bool connected = false;

    /// Solve using the clique algorithm instead?
    bool clique = false;
};

struct CommonSubgraphResult
{
    /// The mapping, empty if none found.
    VertexToVertexMapping mapping;

    /// Total number of nodes processed (recursive calls).
    unsigned long long nodes = 0;

    /// Number of solutions, only if enumerating
    loooong solution_count = 0;

    /// Extra stats, to output
    std::list<std::string> extra_stats;

    /// Did we perform a complete search?
    bool complete = false;
};

auto solve_common_subgraph_problem(
        const InputGraph & first,
        const InputGraph & second,
        const CommonSubgraphParams & params) -> CommonSubgraphResult;

#endif
