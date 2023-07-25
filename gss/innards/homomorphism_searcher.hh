#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_SEARCHER_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_SEARCHER_HH 1

#include <gss/homomorphism.hh>
#include <gss/innards/homomorphism_domain.hh>
#include <gss/innards/homomorphism_model.hh>
#include <gss/innards/homomorphism_traits.hh>
#include <gss/innards/watches.hh>

#include <functional>
#include <random>

namespace gss::innards
{
    enum class SearchResult
    {
        Aborted,
        Unsatisfiable,
        UnsatisfiableAndBackjumpUsingLackey,
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

        auto assignments_as_proof_decisions(const HomomorphismAssignments & assignments) const -> std::vector<std::pair<int, int>>;

        auto solution_in_proof_form(const HomomorphismAssignments & assignments) const -> std::vector<std::pair<NamedVertex, NamedVertex>>;

        template <bool directed_, bool has_edge_labels_, bool induced_, bool verbose_proofs_>
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
            const DuplicateSolutionFilterer &, const std::shared_ptr<Proof> &);

        auto expand_to_full_result(const HomomorphismAssignments & assignments, VertexToVertexMapping & mapping) -> void;

        auto propagate(bool initial, Domains & new_domains, HomomorphismAssignments & assignments, bool propagate_using_lackey) -> bool;

        auto restarting_search(
            HomomorphismAssignments & assignments,
            const Domains & domains,
            unsigned long long & nodes,
            unsigned long long & propagations,
            loooong & solution_count,
            int depth,
            RestartsSchedule & restarts_schedule) -> SearchResult;

        auto save_result(const HomomorphismAssignments & assignments, HomomorphismResult & result) -> void;

        auto set_seed(int n) -> void;

        Watches<HomomorphismAssignment, HomomorphismAssignmentWatchTable> watches;
    };
}

#endif
