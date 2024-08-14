#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_SEARCHER_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_SEARCHER_HH 1

#include <gss/homomorphism.hh>
#include <gss/innards/homomorphism_domain.hh>
#include <gss/innards/homomorphism_model.hh>
#include <gss/innards/homomorphism_traits.hh>
#include <gss/innards/watches.hh>
#include <gss/innards/automorphisms.hh>

#include <functional>
#include <random>

#include "dejavu.h"

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

        // Maybe should be moved into Model?
        std::vector<std::pair<unsigned,unsigned>> useful_target_constraints, useful_pattern_constraints;
        std::vector<int> target_base, pattern_base; 
        std::vector<int> symmetric_value_displacement;
        dejavu::groups::random_schreier t_rschreier{model.target_size * 3}, p_rschreier{model.pattern_size * 3};    // TODO the * 3 is a clunky upper bound in directed cases, needs refining

        auto assignments_as_proof_decisions(const HomomorphismAssignments & assignments) const -> std::vector<std::pair<int, int>>;

        auto solution_in_proof_form(const HomomorphismAssignments & assignments) const -> std::vector<std::pair<NamedVertex, NamedVertex>>;

        template <bool directed_, bool has_edge_labels_, bool induced_, bool verbose_proofs_>
        auto propagate_adjacency_constraints(HomomorphismDomain & d, const HomomorphismAssignment & current_assignment) -> void;

        auto both_in_the_neighbourhood_of_some_vertex(unsigned v, unsigned w) -> bool;

        auto propagate_simple_constraints(Domains & new_domains, const HomomorphismAssignment & current_assignment) -> bool;

        auto propagate_less_thans(Domains & new_domains) -> bool;

        auto propagate_less_thans(Domains & new_domains, const std::vector<std::pair<unsigned int, unsigned int>> & constraints) -> bool;

        auto propagate_less_thans(HomomorphismAssignments assignments, const std::vector<std::pair<unsigned int, unsigned int>> & constraints) -> bool;

        auto propagate_occur_less_thans(const std::optional<HomomorphismAssignment> &, const HomomorphismAssignments &, Domains & new_domains) -> bool;

        auto propagate_dynamic_occur_less_thans(const std::optional<HomomorphismAssignment> &, const HomomorphismAssignments &, Domains & new_domains) -> bool;

        auto make_useful_target_constraints(const std::optional<HomomorphismAssignment> &current_assignment,std::vector<std::pair<unsigned int, unsigned int>> &useful_constraints, std::vector<int> &base) -> bool;
        
        auto make_useful_pattern_constraints(const std::optional<HomomorphismAssignment> &current_assignment,std::vector<std::pair<unsigned int, unsigned int>> &useful_constraints,  std::vector<int> &base) -> bool;

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
            Domains & domains,
            unsigned long long & nodes,
            unsigned long long & propagations,
            loooong & solution_count,
            int depth,
            RestartsSchedule & restarts_schedule,
            std::vector<int> & pattern_orbit_base,
            std::vector<int> & target_orbit_base) -> SearchResult;

        auto save_result(const HomomorphismAssignments & assignments, HomomorphismResult & result) -> void;

        auto set_seed(int n) -> void;

        auto print_pattern_constraints() -> void;
        auto print_target_constraints() -> void;

        int sym_time;

        Watches<HomomorphismAssignment, HomomorphismAssignmentWatchTable> watches;
    };
}

#endif
