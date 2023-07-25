#ifndef GLASGOW_SUBGRAPH_SOLVER_COMMON_SUBGRAPH_HH
#define GLASGOW_SUBGRAPH_SOLVER_COMMON_SUBGRAPH_HH 1

#include <gss/formats/input_graph.hh>
#include <gss/innards/proof-fwd.hh>
#include <gss/loooong.hh>
#include <gss/proof_options.hh>
#include <gss/timeout.hh>
#include <gss/vertex_to_vertex_mapping.hh>

#include <functional>
#include <list>
#include <memory>
#include <string>

namespace gss
{
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
        std::function<auto(const VertexToVertexMapping &)->void> enumerate_callback;

        /// Optional proof options
        std::optional<ProofOptions> proof_options;

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
}

#endif
