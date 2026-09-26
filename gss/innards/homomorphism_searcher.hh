#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_SEARCHER_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_SEARCHER_HH 1

#include <gss/homomorphism.hh>
#include <gss/innards/cost_bound.hh>
#include <gss/innards/filter_activations.hh>
#include <gss/innards/homomorphism_domain.hh>
#include <gss/innards/homomorphism_model.hh>
#include <gss/innards/homomorphism_traits.hh>
#include <gss/innards/watches.hh>

#include <functional>
#include <memory>
#include <optional>
#include <random>

namespace gss::innards
{
    enum class SearchResult
    {
        Aborted,
        Unsatisfiable,
        Satisfiable,
        SatisfiableButKeepGoing,
        Restart
    };

    struct HomomorphismAssignment
    {
        unsigned pattern_vertex;
        unsigned target_vertex;

        auto operator==(const HomomorphismAssignment & other) const -> bool
        {
            return pattern_vertex == other.pattern_vertex && target_vertex == other.target_vertex;
        }

        auto operator!=(const HomomorphismAssignment & other) const -> bool
        {
            return ! (*this == other);
        }
    };

    template <typename EntryType_>
    struct HomomorphismAssignmentWatchTable
    {
        unsigned target_size;
        std::vector<EntryType_> data;

        EntryType_ & operator[](HomomorphismAssignment x)
        {
            return data[target_size * x.pattern_vertex + x.target_vertex];
        }
    };

    struct HomomorphismAssignmentInformation
    {
        HomomorphismAssignment assignment;
        bool is_decision;
        int discrepancy_count;
        int choice_count;
    };

    struct HomomorphismAssignments
    {
        std::vector<HomomorphismAssignmentInformation> values;

        bool contains(const HomomorphismAssignment & assignment) const
        {
            // this should not be a linear scan...
            return values.end() != find_if(values.begin(), values.end(), [&](const auto & a) {
                return a.assignment == assignment;
            });
        }
    };

    using DuplicateSolutionFilterer = const std::function<auto(const HomomorphismAssignments &)->bool>;

    class HomomorphismSearcher
    {
    private:
        using Domains = std::vector<HomomorphismDomain>;

        const HomomorphismModel & model;
        const HomomorphismParams & params;
        const DuplicateSolutionFilterer _duplicate_solution_filterer;

        const std::shared_ptr<Proof> proof;

        std::mt19937 global_rand;

        // Filter-activation accounting (see filter_activations.hh), and whether anything at
        // all wants to know which graph pair removed which value: verbose proof comments do,
        // and so does activation recording. Neither is on in a normal solve, and attributing
        // a removal costs a popcount per graph pair, so this selects the instantiation of
        // propagate_adjacency_constraints() that does the attributing rather than being
        // tested inside its loop.
        const bool _record_filter_activations, _verbose_proof_comments, _track_removals;
        FilterActivations _filter_activations;

        // When minimising cost: the bound, and the best mapping found so far. Every
        // mapping search reaches is cheaper than the last, since the bound refuses
        // anything that is not.
        std::unique_ptr<CostBound> _cost_bound;
        long long _incumbent_cost = CostBound::infinity();
        std::optional<HomomorphismAssignments> _incumbent;

        auto assigned_targets(const HomomorphismAssignments & assignments) const -> std::vector<int>;

        auto assignments_as_proof_decisions(const HomomorphismAssignments & assignments) const -> std::vector<std::pair<int, int>>;

        auto solution_in_proof_form(const HomomorphismAssignments & assignments) const -> std::vector<std::pair<NamedVertex, NamedVertex>>;

        template <bool directed_, bool has_edge_labels_, bool induced_, bool track_removals_>
        auto propagate_adjacency_constraints(HomomorphismDomain & d, const HomomorphismAssignment & current_assignment) -> void;

        auto both_in_the_neighbourhood_of_some_vertex(unsigned v, unsigned w) -> bool;

        auto propagate_simple_constraints(Domains & new_domains, const HomomorphismAssignment & current_assignment) -> bool;

        auto propagate_less_thans(Domains & new_domains) -> bool;

        auto propagate_occur_less_thans(const std::optional<HomomorphismAssignment> &, const HomomorphismAssignments &, Domains & new_domains) -> bool;

        auto find_branch_domain(const Domains & domains) -> const HomomorphismDomain *;

        auto copy_nonfixed_domains_and_make_assignment(
            const Domains & domains,
            unsigned branch_v,
            unsigned f_v) -> Domains;

        auto post_nogood(
            const HomomorphismAssignments & assignments) -> void;

        auto softmax_shuffle(
            std::vector<int> & branch_v,
            unsigned branch_v_end) -> void;

        auto degree_sort(
            std::vector<int> & branch_v,
            unsigned branch_v_end,
            bool reverse) -> void;

    public:
        /**
         * \param cost_data if not null, search for a cheapest mapping (see
         *     HomomorphismParams::minimise_cost). Must outlive the searcher.
         */
        HomomorphismSearcher(const HomomorphismModel & m, const HomomorphismParams & p,
            const DuplicateSolutionFilterer &, const std::shared_ptr<Proof> &,
            Watches<HomomorphismAssignment, HomomorphismAssignmentWatchTable> & watches,
            const CostData * cost_data = nullptr);

        /**
         * When minimising cost, the cheapest mapping found and its cost, if any.
         */
        auto incumbent() const -> const std::optional<HomomorphismAssignments> &;
        auto incumbent_cost() const -> long long;

        auto expand_to_full_result(const HomomorphismAssignments & assignments, VertexToVertexMapping & mapping) -> void;

        auto propagate(bool initial, Domains & new_domains, HomomorphismAssignments & assignments) -> bool;

        auto restarting_search(
            HomomorphismAssignments & assignments,
            const Domains & domains,
            unsigned long long & nodes,
            unsigned long long & propagations,
            loooong & solution_count,
            int depth,
            RestartsSchedule & restarts_schedule) -> SearchResult;

        auto save_result(const HomomorphismAssignments & assignments, HomomorphismResult & result) -> void;

        // Append this searcher's filter-activation counts, if it was recording any.
        auto add_extra_stats(std::list<std::string> & stats) const -> void;

        auto set_seed(int n) -> void;

        // The nogood store is owned by the carried SolveState (sequential path) or
        // by a per-thread vector (threaded path), not by the searcher, so that
        // nogoods can persist across staged search rounds independently of any one
        // searcher's lifetime.
        Watches<HomomorphismAssignment, HomomorphismAssignmentWatchTable> & watches;
    };
}

#endif
