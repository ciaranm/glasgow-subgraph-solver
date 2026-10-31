#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_GSS_INNARDS_SPLIT_PROOF_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_GSS_INNARDS_SPLIT_PROOF_HH 1

#include <gss/innards/proof-fwd.hh>

#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace gss::innards
{
    class HomomorphismProofs;

    // The proof that goes with a counting statement written by a search that splits into
    // components (HomomorphismParams::split_components): from the solver's own encoding of
    // the instance to the statement's root, and then on to the statement's tagged formula,
    // ending EQUIENUMERABLE. See dev_docs/split-counting.md, and counting-statements.md for
    // the format and the tagged formula.
    //
    // The searcher builds the statement as it goes, through the node functions here, each of
    // which defines its node's gate into the core as soon as the node is written. The gates
    // are the proof's own: the tagged formula does not mention them.
    class SplitProof
    {
    private:
        struct Imp;
        std::unique_ptr<Imp> _imp;

        // From the core that the search leaves, the gates' definitions and the root, to
        // exactly the tagged formula.
        auto finish_tagged(int root, long root_line) -> void;

    public:
        SplitProof(const std::shared_ptr<Proof> &, HomomorphismProofs &, int pattern_size, int target_size,
            const std::string & statement_file);

        // Removes the statement file, unless finish() wrote it in full.
        ~SplitProof();

        SplitProof(const SplitProof &) = delete;
        auto operator=(const SplitProof &) -> SplitProof & = delete;

        // Copy each assignment variable to the statement's own, and move the preserved set
        // across. Call once, at level 0, before searching.
        auto begin() -> void;

        // Statement nodes. A conjunction or disjunction of one node is that node. The
        // children of a disjunction are told apart by one pattern vertex: each child puts it
        // at its value, distinct from the others', except that at most one child may have
        // the value -1, meaning that it puts the vertex at none of the others' values. The
        // proof of the statement's count relies on this.
        auto cube(std::vector<std::pair<int, int>> assignment) -> int;
        auto conjunction(std::vector<int> children) -> int;
        auto disjunction(const std::vector<int> & children, int vertex, const std::vector<int> & values) -> int;

        auto gate(int node) const -> std::string;
        auto statement_variable(int v, int t) const -> std::string;

        // A fresh literal, defined at level 0, for none of these (pattern vertex, target
        // vertex) assignments being made.
        auto assumption(const std::vector<std::pair<int, int>> & excluded) -> std::string;

        // "rup <literals> >= 1", at the current level or at level 0, returning its line.
        auto rup(const std::vector<std::string> & literals, bool at_level_zero) -> long;

        // For each (node, pattern vertex, target vertex), derive that the node's gate puts
        // the vertex there, which must hold of every cube below the node.
        auto derive_signatures(const std::vector<std::tuple<int, int, int>> & lemmas) -> void;

        // Delete these derived lines.
        auto forget(const std::vector<long> & lines) -> void;

        // The search found solutions, and the proof has derived that the root holds. Turn
        // the core into exactly the statement's tagged formula, and conclude that it is
        // equienumerable with the input.
        auto finish(int root) -> void;
    };
}

#endif
