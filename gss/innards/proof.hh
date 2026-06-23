#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_PROOF_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_PROOF_HH 1

#include <gss/innards/proof-fwd.hh>
#include <gss/loooong.hh>
#include <gss/proof_options.hh>

#include <exception>
#include <functional>
#include <iosfwd>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace gss::innards
{
    class ProofError : public std::exception
    {
    private:
        std::string _message;

    public:
        ProofError(const std::string & message) noexcept;

        virtual auto what() const noexcept -> const char *;
    };

    using NamedVertex = std::pair<int, std::string>;

    // The proof state shared across the homomorphism adjacency / supplemental-graph
    // derivations: for each (graph g, pattern p, pattern q, target t) the adjacency
    // constraint's label, its numeric line id (for deleting a transiently-derived one),
    // and its permitted target set (needed to cancel a self-loop term). These are keyed
    // entirely by index, not name. Grouping them lets the derivations that read and write
    // them migrate out of Proof into the solver-proofs middle layer one at a time, each
    // sharing this one object rather than a Proof-private map. (See
    // dev_docs/preprocessor-refactor.md, Phase 3 Option 2.)
    struct AdjacencyProofLines
    {
        std::map<std::tuple<long, long, long, long>, std::string> labels;
        std::map<std::tuple<long, long, long, long>, long> ids;
        std::map<std::tuple<long, long, long, long>, std::vector<long>> permitted;
    };

    class Proof
    {
    private:
        struct Imp;
        std::unique_ptr<Imp> _imp;

        auto need_elimination(int p, int t) -> void;

    public:
        explicit Proof(const ProofOptions &);
        Proof(Proof &&);
        ~Proof();
        auto operator=(Proof &&) -> Proof &;

        Proof(const Proof &) = delete;
        auto operator=(const Proof &) -> Proof & = delete;

        auto super_extra_verbose() const -> bool;

        // model-writing functions
        auto create_cp_variable(int pattern_vertex, int target_size,
            const std::function<auto(int)->std::string> & pattern_name,
            const std::function<auto(int)->std::string> & target_name) -> void;

        auto create_injectivity_constraints(int pattern_size, int target_size,
            const std::function<auto(int)->std::string> & target_name) -> void;

        // Local-injectivity analogue: for each pattern vertex v and target t, at most one
        // neighbour of v maps to t (phi restricted to N(v) is injective).
        auto create_locally_injective_constraints(int pattern_size, int target_size,
            const std::function<auto(int, int)->bool> & adjacent,
            const std::function<auto(int)->std::string> & pattern_name,
            const std::function<auto(int)->std::string> & target_name) -> void;

        auto create_forbidden_assignment_constraint(int p, int t) -> void;
        auto start_adjacency_constraints_for(int p, int t) -> void;
        auto create_adjacency_constraint(const NamedVertex & p, const NamedVertex & q, const NamedVertex & t,
            const std::vector<int> & u, bool induced) -> void;

        // Generic low-level proof-emission primitives, so the homomorphism-specific
        // derivations can live in the solver-proofs middle layer and emit through here
        // rather than reaching into Proof's innards. emit_proof_line writes a constraint
        // line (appending the newline), bumps the proof-line counter, and returns the new
        // line's number; emit_proof_directive writes a non-constraint line (a comment, a
        // setlvl/wiplvl, a `del id`) with no counter change; current_proof_line is the
        // most recent line's number, for building a citation. variable_name is the proof
        // variable for assigning pattern vertex p to target vertex t.
        auto emit_proof_line(const std::string & line) -> long;
        auto emit_proof_directive(const std::string & line) -> void;
        [[nodiscard]] auto current_proof_line() const -> long;
        [[nodiscard]] auto variable_name(int p, int t) const -> const std::string &;
        // True when create_cp_variable recorded a variable for (p, t).  Used by
        // the POL loop-cancellation to skip pattern nodes whose domain was pruned
        // before injectivity constraints were written.
        [[nodiscard]] auto has_variable_mapping(int p, int t) const -> bool;
        [[nodiscard]] auto is_locally_injective() const -> bool;

        // The OPB-model analogues: emit_model_constraint writes a constraint into the model
        // (bumping the constraint count); emit_model_comment writes a `*` comment line.
        auto emit_model_constraint(const std::string & line) -> void;
        auto emit_model_comment(const std::string & line) -> void;

        // Read accessors for the model constraint labels the homomorphism derivations cite:
        // the injectivity constraint on target t, its local-injectivity analogue on (pattern
        // p, target t), and the at-most-one-value constraint on pattern vertex p. Plus the
        // generic dedup cache (keyed by a constraint's text) the supplemental derivations use
        // to reuse an identical line's label instead of re-deriving it.
        [[nodiscard]] auto injectivity_label(int t) const -> const std::string &;
        [[nodiscard]] auto locally_injective_label(int p, int t) const -> const std::string &;
        [[nodiscard]] auto at_most_one_value_label(int p) const -> const std::string &;
        [[nodiscard]] auto cached_proof_line(const std::string & key) const -> std::optional<std::string>;
        auto cache_proof_line(const std::string & key, const std::string & label) -> void;

        // Declare a projected (preserved) set, listing exactly the assignment
        // variables, so the proof's solution count is in terms of the high-level
        // mapping rather than any auxiliary encoding variables. Must be called
        // after all create_cp_variable calls but before finalise_model.
        auto emit_preserved_assignment_variables() -> void;

        auto finalise_model() -> void;

        // when we're done
        auto finish_unsat_proof() -> void;
        auto finish_sat_proof() -> void;
        auto finish_unknown_proof() -> void;
        auto finish_optimisation_proof(int size) -> void;

        // Conclude a counting / enumeration proof: ENUMERATION_COMPLETE if the
        // whole search space was exhausted, otherwise ENUMERATION_PARTIAL.
        auto finish_enumeration_proof(const loooong & number_of_solutions, bool complete) -> void;

        // top of search failures
        auto failure_due_to_pattern_bigger_than_target() -> void;

        // domain initialisation
        auto incompatible_by_degrees(
            int g,
            const NamedVertex & p,
            const std::vector<int> & n_p,
            const NamedVertex & t,
            const std::vector<int> & n_t) -> void;

        auto incompatible_by_nds(
            int g,
            const NamedVertex & p,
            const NamedVertex & t,
            const std::vector<int> & p_subsequence,
            const std::vector<int> & t_subsequence,
            const std::vector<int> & t_remaining) -> void;

        auto incompatible_by_loops(
            const NamedVertex & p,
            const NamedVertex & t) -> void;

        auto initial_domain_is_empty(int p, const std::string & where) -> void;

        // The shared adjacency-line proof state, so the derivations that build and consume it
        // can move into the solver-proofs middle layer while it still lives here.
        [[nodiscard]] auto adjacency_proof_lines() -> AdjacencyProofLines &;

        // new constraints
        auto emit_hall_set_or_violator(const std::vector<NamedVertex> & lhs, const std::vector<NamedVertex> & rhs) -> void;

        // branch logging
        auto root_propagation_failed() -> void;
        auto guessing(int depth, const NamedVertex & branch_v, const NamedVertex & val) -> void;
        auto propagation_failure(const std::vector<std::pair<int, int>> & decisions, const NamedVertex & branch_v, const NamedVertex & val) -> void;
        auto incorrect_guess(const std::vector<std::pair<int, int>> & decisions, bool was_failure) -> void;
        auto out_of_guesses(const std::vector<std::pair<int, int>> & decisions) -> void;
        auto unit_propagating(const NamedVertex & var, const NamedVertex & val) -> void;

        // proof levels
        auto start_level(int level) -> void;
        auto back_up_to_level(int level) -> void;
        auto forget_level(int level) -> void;
        auto back_up_to_top() -> void;
        auto post_restart_nogood(const std::vector<std::pair<int, int>> & decisions) -> void;

        // cliques
        auto create_binary_variable(int vertex,
            const std::function<auto(int)->std::string> & name) -> void;
        auto create_objective(int n, std::optional<int> d) -> void;
        auto create_non_edge_constraint(const NamedVertex & p, const NamedVertex & q) -> void;
        auto backtrack_from_binary_variables(const std::vector<int> &) -> void;
        auto colour_bound(const std::vector<std::vector<int>> &) -> void;

        // clique for hom
        auto prepare_hom_clique_proof(const NamedVertex & p,
            const NamedVertex & t,
            unsigned size) -> void;

        auto start_hom_clique_proof(const NamedVertex & p,
            std::vector<NamedVertex> && p_clique_neighbourhood,
            const NamedVertex & t,
            std::map<int, NamedVertex> && t_clique_neighbourhood) -> void;

        auto finish_hom_clique_proof(const NamedVertex & p, const NamedVertex & t, unsigned size) -> void;

        auto add_hom_clique_non_edge(
            const NamedVertex & filter_p,
            const NamedVertex & filter_t,
            const std::vector<NamedVertex> & p_clique,
            const NamedVertex & t,
            const NamedVertex & u) -> void;

        // common subgraphs
        auto create_null_decision_bound(int p, int t, std::optional<int> d) -> void;
        auto mcs_bound(
            const std::vector<std::pair<std::set<int>, std::set<int>>> & partitions) -> void;
        auto create_connected_constraints(int p, int t, const std::function<auto(int, int)->bool> & adj) -> void;

        // common subgraph to clique
        auto has_clique_model() const -> bool;
        auto create_clique_encoding(const std::vector<std::pair<int, int>> &, const std::vector<std::pair<int, int>> &) -> void;
        auto create_clique_nonedge(int v, int w) -> void;
        auto not_connected_in_underlying_graph(const std::vector<int> &, int) -> void;

        // enumeration
        auto post_solution(const std::vector<std::pair<NamedVertex, NamedVertex>> & decisions) -> void;
        auto post_solution(const std::vector<int> & solution) -> void;

        // optimisation
        auto new_incumbent(const std::vector<std::pair<int, bool>> & solution) -> void;
        auto new_incumbent(const std::vector<std::tuple<NamedVertex, NamedVertex, bool>> & solution) -> void;

        // super extra verbose
        auto show_domains(const std::string & where, const std::vector<std::pair<NamedVertex, std::vector<NamedVertex>>> & domains) -> void;
        auto propagated(const NamedVertex & p, const NamedVertex & t, int g, int n_values, const NamedVertex & q) -> void;
    };
}

#endif
