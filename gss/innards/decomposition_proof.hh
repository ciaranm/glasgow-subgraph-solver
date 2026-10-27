#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_GSS_INNARDS_DECOMPOSITION_PROOF_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_GSS_INNARDS_DECOMPOSITION_PROOF_HH 1

#include <gss/counting_statement.hh>
#include <gss/formats/input_graph.hh>
#include <gss/homomorphism.hh>

#include <map>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace gss::innards
{
    // The proof that goes with a counting statement written by --decomposition: from the
    // solver's own encoding of the instance to the statement's root. This is phases 4 and 5
    // of dev_docs/decomposition-counting.md, and the steps it takes are described there and
    // in dev_docs/counting-statements.md ("From an input to the counting formula").
    //
    // It writes the model for the whole instance, then copies each assignment variable to
    // the statement's own variable and moves the preserved set across. Every node of the
    // statement gets its gate defined, into the core, as soon as it is written (see
    // define_gate, which is meant to be a CountingStatementWriter observer). Each component
    // is searched under proof, logging into this proof under a guard that folds in each
    // solution's cube as it is found, so that the search concludes that the component's
    // leaf holds.
    //
    // With hubs, the caller walks the tree of guesses, and each component search runs under
    // an assumption: that its vertices are off the hubs, and that its neighbours on hubs
    // are where the guess puts them. Everything this needs that is not in the statement is
    // scaffolding: derived, never in the core, and deleted once the root is derived.
    class DecompositionProof
    {
    private:
        struct Imp;
        std::unique_ptr<Imp> _imp;

    public:
        DecompositionProof(const ProofOptions &, const InputGraph & pattern, const InputGraph & target, const std::vector<int> & hubs);
        ~DecompositionProof();

        DecompositionProof(const DecompositionProof &) = delete;
        auto operator=(const DecompositionProof &) -> DecompositionProof & = delete;

        // Define node i's gate, both directions, and put the definition in the core.
        auto define_gate(int i, const CountingNode & node) -> void;

        struct SearchedComponent
        {
            // each embedding, in the component's own numbering
            std::vector<std::vector<int>> embeddings;
            std::vector<int> cubes;
            bool complete = false;
            // the search's conclusion: the last fold, or the assumption failing if there
            // were no embeddings, and the folds' definitions
            std::optional<long> conclusion;
            std::vector<long> fold_lines;
        };

        // Search one connected component of the pattern under proof, under the assumption
        // if there is one, writing a cube to the statement for each embedding as it is
        // found, or reusing the cube for that assignment if it has already been written.
        // The component's vertex names must be the pattern's names for the same vertices,
        // which are vertices[0], vertices[1], ... If there are no embeddings, the proof has
        // derived that the assumption fails, or a contradiction if there is none.
        auto search_component(const InputGraph & component, const std::vector<int> & vertices,
            HomomorphismParams params, CountingStatementWriter & writer,
            const std::optional<std::string> & assumption) -> SearchedComponent;

        // After a search that found embeddings, derive that the leaf, which must contain
        // every cube the search found, holds under the assumption, and forget the folds.
        auto conclude_component(const SearchedComponent &, int leaf, const std::optional<std::string> & assumption) -> void;

        // Derive a literal, by unit propagation, at level 0, returning its line.
        auto derive(const std::string & literal) -> long;

        // Derive a clause, by unit propagation, at level 0, as scaffolding, returning its
        // line. With no literals, this is a contradiction.
        auto derive_clause(const std::vector<std::string> & literals) -> long;

        // Delete these scaffolding lines now.
        auto forget(const std::vector<long> & lines) -> void;

        // Delete all of the scaffolding, once the root has been derived.
        auto forget_scaffolding() -> void;

        // The literal for pattern vertex v going to target vertex t.
        auto assignment(int v, int t) const -> std::string;

        // The literal for pattern vertex v going to no hub, defined when first asked for.
        auto off_hubs(int v) -> std::string;

        // A literal for all of these literals holding, defined when first asked for.
        auto conjunction(std::vector<std::string> literals) -> std::string;

        // Pattern vertex v cannot go to non-hub t if its neighbours w are off the hubs,
        // because they would need more distinct non-hub neighbours of t than t has.
        auto restrict_by_degree(int v, int t, const std::vector<int> & w) -> void;

        // Pattern vertex v cannot go to target vertex t, by unit propagation once it is
        // assumed to, for instance because v has a loop and t does not. (The other ways a
        // vertex cannot go on a hub follow from the decisions already made, without this.)
        auto rule_out(int v, int t) -> void;

        // Pattern vertex v cannot go to hub h, because it has more neighbours than h does.
        auto hub_too_small(int v, int h) -> void;

        // A sum that says these pattern vertices cannot all be off the hubs, because there
        // are more of them than there are non-hubs. The caller derives what it needs from
        // this by unit propagation.
        auto too_many_off_hubs(const std::vector<int> & vertices) -> void;

        // Fold in a guess's gate: the new fold literal holds if and only if the previous
        // one does, or the gate does.
        auto fold_guess(int gate) -> std::string;

        // Every cube of a glued leaf sends two pattern vertices to one target vertex, so
        // none of them, and hence not the leaf, can hold.
        auto refute_glued_leaf(int leaf, const std::vector<int> & cubes) -> void;

        // The instance has no solutions: the pattern has more vertices than the target.
        auto finish_pattern_bigger_than_target() -> void;

        // The instance has no solutions, and the proof has derived a contradiction.
        auto finish_no_solutions() -> void;

        // What finish needs to know about the statement, once its root has been derived.
        struct Structure
        {
            int root;
            long root_line;

            // each component of what a guess leaves off the hubs: its vertices, its leaf,
            // and its cubes, each with its embedding as the image of each vertex in turn
            struct Component
            {
                std::vector<int> vertices;
                int leaf;
                std::vector<int> cubes;
                std::vector<std::vector<int>> embeddings;
            };

            // each guess under the root, or just the root if it is the only one: its gate,
            // the or of the overlaps of what it leaves off the hubs, if there are any, the
            // components of that, and each overlap by its partition of those vertices
            // (sorted blocks, sorted, singletons included)
            struct Guess
            {
                int gate;
                std::optional<int> overlaps;
                std::vector<Component> components;
                std::map<std::vector<std::vector<int>>, int> overlaps_by_partition;
            };
            std::vector<Guess> guesses;
        };

        // Delete everything but the counting formula from the core, which shows that the
        // root implies the encoding, and conclude that the core is equienumerable with
        // the translator's counting formula (phase 5 of decomposition-counting.md).
        auto finish(const Structure &) -> void;
    };
}

#endif
