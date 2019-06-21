/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_HOMOMORPHISM_HH
#define GLASGOW_SUBGRAPH_SOLVER_HOMOMORPHISM_HH 1

#include "formats/input_graph.hh"
#include "lackey.hh"
#include "restarts.hh"
#include "timeout.hh"
#include "value_ordering.hh"
#include "vertex_to_vertex_mapping.hh"

#include <functional>
#include <list>
#include <memory>
#include <string>

enum class Injectivity
{
    Injective,
    LocallyInjective,
    NonInjective
};

struct HomomorphismParams
{
    /// Timeout handler
    std::shared_ptr<Timeout> timeout;

    /// The start time of the algorithm.
    std::chrono::time_point<std::chrono::steady_clock> start_time;

    /// Induced?
    bool induced = false;

    /// Noninjective?
    Injectivity injectivity = Injectivity::Injective;

    /// Enumerate?
    bool count_solutions = false;

    /// Find a minimal unsat pattern?
    bool minimal_unsat_pattern = false;

    /// Print solutions, for enumerating
    std::function<auto (const VertexToVertexMapping &) -> void> enumerate_callback;

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

    /// Are we allowed to do clique detection?
    bool clique_detection = true;

    /// Are we allowed to remove isolated vertices?
    bool remove_isolated_vertices = true;

    /// Use common neighbour shape filtering?
    bool common_neighbour_shapes = false;

    /// Less pattern constraints
    std::list<std::pair<std::string, std::string> > pattern_less_constraints;

    /// Optional lackey, for external side constraints
    std::unique_ptr<Lackey> lackey;
};

struct HomomorphismResult
{
    /// The mapping, empty if none found.
    VertexToVertexMapping mapping;

    /// Vertices in a minimal unsat pattern, if requested.
    std::list<int> minimal_unsat_pattern;

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
