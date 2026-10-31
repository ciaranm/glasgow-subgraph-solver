#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_SEARCHER_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_SEARCHER_HH 1

#include <gss/homomorphism.hh>
#include <gss/innards/filter_activations.hh>
#include <gss/innards/homomorphism_domain.hh>
#include <gss/innards/homomorphism_model.hh>
#include <gss/innards/homomorphism_traits.hh>
#include <gss/innards/watches.hh>

#include <cstddef>
#include <functional>
#include <map>
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

    struct ComponentTally;
    struct SplitCaches;
    class SplitProof;

    // With a counting statement (see SplitProof), what a search node hands back to its
    // parent: for each signature of the solutions below it (empty when not inside a split's
    // component), the statement node they lie in, over the vertices unassigned at the node,
    // or -1 if there are none of those and the node is a solution.
    using StatementEntries = std::vector<std::pair<std::vector<int>, int>>;

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

        // While one component of a split is being searched, its solutions go into this tally
        // rather than into the solution count (see split_into_components). Null otherwise.
        ComponentTally * _tally = nullptr;
        unsigned long long _splits = 0, _nested_splits = 0;
        std::unique_ptr<SplitCaches> _split_caches;

        // With a counting statement: what the last search node returned, and where the
        // decisions of the component search under way begin in the assignments, so that its
        // lines mention only its own decisions and its assumption.
        SplitProof * _split_proof = nullptr;
        StatementEntries _entries;
        std::size_t _decision_base = 0;

        // "rup", the negated decisions, the assumption's negation if there is one, these
        // literals, at the current level.
        auto statement_line(const HomomorphismAssignments & assignments, const std::vector<std::string> & literals) -> long;

        // After a child of a search node has returned, with the child's assignments still
        // made: add its statement nodes to the node's, by signature, and derive that the
        // child's decision leads to one of them.
        auto add_contributions(const HomomorphismAssignments & assignments, std::size_t assignments_size,
            std::map<std::vector<int>, std::vector<std::pair<int, int>>> & contributions) -> void;

        auto split_into_components(
            HomomorphismAssignments & assignments,
            const Domains & domains,
            unsigned long long & nodes,
            unsigned long long & propagations,
            loooong & solution_count,
            int depth,
            RestartsSchedule & restarts_schedule) -> std::optional<SearchResult>;

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
        HomomorphismSearcher(const HomomorphismModel & m, const HomomorphismParams & p,
            const DuplicateSolutionFilterer &, const std::shared_ptr<Proof> &,
            Watches<HomomorphismAssignment, HomomorphismAssignmentWatchTable> & watches);
        ~HomomorphismSearcher();

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

        // Build a counting statement as the search goes, and prove it, instead of logging each
        // solution. After a complete search, statement_entries() is the root's.
        auto set_split_proof(SplitProof *) -> void;
        auto statement_entries() const -> const StatementEntries &;

        // The nogood store is owned by the carried SolveState (sequential path) or
        // by a per-thread vector (threaded path), not by the searcher, so that
        // nogoods can persist across staged search rounds independently of any one
        // searcher's lifetime.
        Watches<HomomorphismAssignment, HomomorphismAssignmentWatchTable> & watches;
    };
}

#endif
