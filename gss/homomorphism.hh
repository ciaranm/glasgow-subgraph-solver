#ifndef GLASGOW_SUBGRAPH_SOLVER_HOMOMORPHISM_HH
#define GLASGOW_SUBGRAPH_SOLVER_HOMOMORPHISM_HH 1

#include <gss/formats/input_graph.hh>
#include <gss/innards/lackey.hh>
#include <gss/innards/proof-fwd.hh>
#include <gss/loooong.hh>
#include <gss/proof_options.hh>
#include <gss/restarts.hh>
#include <gss/timeout.hh>
#include <gss/value_ordering.hh>
#include <gss/vertex_to_vertex_mapping.hh>

#include <functional>
#include <list>
#include <memory>
#include <optional>
#include <string>

namespace gss
{
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

        /// Noninjective?
        Injectivity injectivity = Injectivity::Injective;

        /// Enumerate?
        bool count_solutions = false;

        /// Print solutions, for enumerating
        std::function<auto(const VertexToVertexMapping &)->bool> enumerate_callback;

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

        /// Any extra shapes to use, with injectivity and count
        std::list<std::tuple<std::unique_ptr<InputGraph>, bool, int>> extra_shapes;

        /// Are we allowed to do clique size constraints?
        bool clique_size_constraints = false;

        /// If we do clique constraints, do we do them on supplemental graphs too?
        bool clique_size_constraints_on_supplementals = false;

        /// Disable neighbourhood degree sequence processing?
        bool no_nds = false;

        /// Less pattern constraints
        std::list<std::pair<std::string, std::string>> pattern_less_constraints;

        /// Occurs less target constraints
        std::list<std::pair<std::string, std::string>> target_occur_less_constraints;

        /// Optional lackey, for external side constraints
        std::unique_ptr<innards::Lackey> lackey;

        /// Send partial solutions to the lackey?
        bool send_partials_to_lackey = false;

        /// Propagate using the lackey?
        PropagateUsingLackey propagate_using_lackey = PropagateUsingLackey::Never;

        /// Optional proof options
        std::optional<ProofOptions> proof_options;
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

        /// Did we perform a complete search?
        bool complete = false;
    };

    auto solve_homomorphism_problem(
        const InputGraph & pattern,
        const InputGraph & target,
        const HomomorphismParams & params) -> HomomorphismResult;
}

#endif
