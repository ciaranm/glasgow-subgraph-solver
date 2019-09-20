/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_PROOF_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_PROOF_HH 1

#include "proof-fwd.hh"
#include "vertex_to_vertex_mapping.hh"

#include <exception>
#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

class ProofError : public std::exception
{
    private:
        std::string _message;

    public:
        ProofError(const std::string & message) noexcept;

        virtual auto what() const noexcept -> const char *;
};

class Proof
{
    private:
        struct Imp;
        std::unique_ptr<Imp> _imp;

    public:
        Proof(const std::string & opb_file, const std::string & log_file, bool levels, bool solutions);
        Proof(Proof &&);
        ~Proof();
        auto operator= (Proof &&) -> Proof &;

        Proof(const Proof &) = delete;
        auto operator= (const Proof &) -> Proof & = delete;

        // model-writing functions
        auto create_cp_variable(int pattern_vertex, int target_size) -> void;
        auto create_injectivity_constraints(int pattern_size, int target_size) -> void;
        auto create_forbidden_assignment_constraint(int p, int t) -> void;
        auto create_adjacency_constraint(int p, int q, int t, const std::vector<int> & u) -> void;
        auto finalise_model() -> void;

        // when we're done
        auto finish_unsat_proof() -> void;

        // top of search failures
        auto failure_due_to_pattern_bigger_than_target() -> void;

        // domain initialisation
        auto incompatible_by_degrees(int g, int p, const std::vector<int> & n_p, int t, const std::vector<int> & n_t) -> void;
        auto incompatible_by_nds(int g, int p, int t) -> void;
        auto initial_domain_is_empty(int p) -> void;

        // new constraints
        auto emit_hall_set_or_violator(const std::vector<int> & lhs, const std::vector<int> & rhs) -> void;

        // branch logging
        auto root_propagation_failed() -> void;
        auto guessing(int depth, int branch_v, int val) -> void;
        auto propagation_failure(const std::vector<std::pair<int, int> > & decisions, int branch_v, int val) -> void;
        auto incorrect_guess(const std::vector<std::pair<int, int> > & decisions) -> void;
        auto out_of_guesses(const std::vector<std::pair<int, int> > & decisions) -> void;
        auto unit_propagating(int var, int val) -> void;

        // proof levels
        auto start_level(int level) -> void;
        auto back_up_to_level(int level) -> void;
        auto back_up_to_top() -> void;
        auto post_restart_nogood(const std::vector<std::pair<int, int> > & decisions) -> void;

        // enumeration
        auto post_solution(const VertexToVertexMapping &, const std::vector<std::pair<int, int> > & decisions) -> void;
};

#endif
