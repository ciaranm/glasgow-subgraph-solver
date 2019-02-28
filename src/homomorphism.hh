/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_HOMOMORPHISM_HH
#define GLASGOW_SUBGRAPH_SOLVER_HOMOMORPHISM_HH 1

#include "formats/input_graph.hh"
#include "restarts.hh"
#include "timeout.hh"
#include "value_ordering.hh"

struct HomomorphismParams
{
    /// Timeout handler
    std::unique_ptr<Timeout> timeout;

    /// The start time of the algorithm.
    std::chrono::time_point<std::chrono::steady_clock> start_time;

    /// Induced?
    bool induced = false;

    /// Noninjective?
    bool noninjective = false;

    /// Enumerate?
    bool enumerate = false;

    /// Which value-ordering heuristic?
    ValueOrdering value_ordering_heuristic = ValueOrdering::Biased;

    /// Restarts schedule
    std::unique_ptr<RestartsSchedule> restarts_schedule;

    /// Largest size of nogood to store (0 disables nogoods)
    unsigned nogood_size_limit = std::numeric_limits<unsigned>::max();

    /// How many threads to use (1 for sequential, 0 to auto-detect). Must be
    /// used in conjunction with restarts.
    unsigned n_threads = 1;

    /// Do one restart before launching remaining threads?
    bool delay_thread_creation = false;

    /// Trigger restarts using the first thread?
    bool triggered_restarts = false;
};

struct HomomorphismResult
{
    /// The mapping, empty if none found.
    std::map<int, int> mapping;

    /// Total number of nodes processed (recursive calls).
    unsigned long long nodes = 0;

    /// Number of times propagate called.
    unsigned long long propagations = 0;

    /// Extra stats, to output
    std::list<std::string> extra_stats;

    /// Number of solutions, only if enumerating
    unsigned long long solution_count = 0;

    /// Did we perform a complete search?
    bool complete = false;
};

auto solve_homomorphism_problem(const std::pair<InputGraph, InputGraph> & graphs, const HomomorphismParams & params) -> HomomorphismResult;

#endif
