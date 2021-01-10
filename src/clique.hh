/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_CLIQUE_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_CLIQUE_HH 1

#include "formats/input_graph.hh"
#include "restarts.hh"
#include "timeout.hh"
#include "proof-fwd.hh"
#include "svo_bitset.hh"

#include <chrono>
#include <functional>
#include <list>
#include <memory>
#include <optional>
#include <set>

enum class ColourClassOrder
{
    ColourOrder,
    SingletonsFirst,
    Sorted
};

struct CliqueParams
{
    /// Timeout handler
    std::shared_ptr<Timeout> timeout;

    /// The start time of the algorithm.
    std::chrono::time_point<std::chrono::steady_clock> start_time;

    /// Decide instead of maximise?
    std::optional<unsigned> decide;

    /// Can stop after finding this size
    std::optional<unsigned> stop_after_finding;

    /// Restarts schedule
    std::unique_ptr<RestartsSchedule> restarts_schedule;

    /// Largest size of nogood to store (0 disables nogoods)
    unsigned nogood_size_limit = std::numeric_limits<unsigned>::max();

    /// Which colour order to use?
    ColourClassOrder colour_class_order = ColourClassOrder::SingletonsFirst;

    /// Colour in input order, rather than degree order
    bool input_order = false;

    /// For use by the maximum common connected subgraph reduction
    std::function<auto (int, const std::function<auto (int) -> int> &) -> SVOBitset> connected;

    /// Optional proof handler
    std::shared_ptr<Proof> proof;

    /// If logging proofs, only log the bound (for use by homomorphism solver for clique filtering)
    bool proof_is_for_hom = false;
};

struct CliqueResult
{
    /// The vertices in the clique, empty if none found.
    std::set<int> clique;

    /// Total number of nodes processed (recursive calls).
    unsigned long long nodes = 0, find_nodes = 0, prove_nodes = 0;

    /// Extra stats, to output
    std::list<std::string> extra_stats;

    /// Did we perform a complete search?
    bool complete = false;
};

auto solve_clique_problem(const InputGraph & graph, const CliqueParams & params) -> CliqueResult;

#endif
