/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_HOMOMORPHISM_HH
#define GLASGOW_SUBGRAPH_SOLVER_HOMOMORPHISM_HH 1

#include "formats/input_graph.hh"
#include "lackey.hh"
#include "loooong.hh"
#include "restarts.hh"
#include "timeout.hh"
#include "value_ordering.hh"
#include "vertex_to_vertex_mapping.hh"
#include "proof-fwd.hh"

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

enum class PropagateUsingLackey
{
    Never,
    Root,
    Always,
    RootAndBackjump
};

struct HomomorphismParams
{
    /// Timeout handler
    std::shared_ptr<Timeout> timeout;

    /// The start time of the algorithm.
    std::chrono::time_point<std::chrono::steady_clock> start_time;

    /// Induced?
    bool induced = false;

    /// Is a bigraph?
    bool bigraph = false;

    /// Are we checking for bigraph equality or just matching?
    bool equality_check = false;

    /// If it is a bigraph, use nogood trickery for projection?
    bool use_bigraph_projection_nogoods = true;

    /// Do the input bigraphs use directed hyperedges?
    bool directed = false;

    /// Noninjective?
    Injectivity injectivity = Injectivity::Injective;

    /// Enumerate?
    bool count_solutions = false;

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

    /// Use distance 3 filtering?
    bool distance3 = false;

    /// Use k4 filtering?
    bool k4 = false;

    /// Disable all supplemental graphs?
    bool no_supplementals = false;

    /// How many exact path graphs do we have, if we have any?
    int number_of_exact_path_graphs = 4;

    /// Disable neighbourhood degree sequence processing?
    bool no_nds = false;

    /// Less pattern constraints
    std::list<std::pair<std::string, std::string> > pattern_less_constraints;

    /// Optional lackey, for external side constraints
    std::unique_ptr<Lackey> lackey;

    /// Send partial solutions to the lackey?
    bool send_partials_to_lackey = false;

    /// Propagate using the lackey?
    PropagateUsingLackey propagate_using_lackey = PropagateUsingLackey::Never;

    /// Optional proof handler
    std::unique_ptr<Proof> proof;
};

struct HomomorphismResult
{
    /// The mapping, empty if none found.
    VertexToVertexMapping mapping;

    /// Total number of nodes processed (recursive calls).
    unsigned long long nodes = 0;

    /// Number of times propagate called.
    unsigned long long propagations = 0;

    /// Extra stats, to output
    std::list<std::string> extra_stats;

    /// Number of solutions, only if enumerating
    loooong solution_count = 0;

    /// Number of solutions that were rejected by the final bigraph check.
    loooong rejected_solution_count = 0;

    /// Did we perform a complete search?
    bool complete = false;
};

auto solve_homomorphism_problem(
        const InputGraph & pattern,
        const InputGraph & target,
        const HomomorphismParams & params) -> HomomorphismResult;

#endif
