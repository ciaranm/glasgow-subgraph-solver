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
    // solver's own encoding of the instance to the statement's root. This is phase 4 of
    // dev_docs/decomposition-counting.md, and the steps it takes are described there and
    // in dev_docs/counting-statements.md ("From an input to the counting formula").
    //
    // It writes the model for the whole instance, then copies each assignment variable to
    // the statement's own variable and moves the preserved set across. Every node of the
    // statement gets its gate defined, into the core, as soon as it is written (see
    // define_gate, which is meant to be a CountingStatementWriter observer). Each component
    // of the pattern is searched under proof, logging into this proof under a guard that
    // folds in each solution's cube as it is found, so that the search concludes that the
    // component's leaf holds.
    class DecompositionProof
    {
    private:
        struct Imp;
        std::unique_ptr<Imp> _imp;

    public:
        DecompositionProof(const ProofOptions &, const InputGraph & pattern, const InputGraph & target);
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
            std::optional<int> leaf;
            bool complete = false;
        };

        // Search one connected component of the pattern under proof, writing a cube to the
        // statement for each embedding as it is found and then the leaf over them, and
        // derive the leaf's gate as a unit. The component's vertex names must be the
        // pattern's names for the same vertices, which are vertices[0], vertices[1], ...
        // If there are no embeddings, nothing but the search is written, and the proof has
        // derived a contradiction.
        auto search_component(const InputGraph & component, const std::vector<int> & vertices,
            HomomorphismParams params, CountingStatementWriter & writer) -> SearchedComponent;

        // Derive a literal, by unit propagation, at level 0, returning its line.
        auto derive(const std::string & literal) -> long;

        // Every cube of a glued leaf sends two pattern vertices to one target vertex, so
        // none of them, and hence not the leaf, can hold.
        auto refute_glued_leaf(int leaf, const std::vector<int> & cubes) -> void;

        // The instance has no solutions: the pattern has more vertices than the target.
        auto finish_pattern_bigger_than_target() -> void;

        // The instance has no solutions, because a component's search derived a
        // contradiction.
        auto finish_no_solutions() -> void;

        // What finish needs to know about the statement, once its root has been derived.
        struct Structure
        {
            int root;
            long root_line;
            // the root's product node, and the or of its overlaps, if it has any
            int product;
            std::optional<int> overlaps;

            // each component of the pattern: its vertices, its leaf, and its cubes, each
            // with its embedding as the image of each vertex in turn
            struct Component
            {
                std::vector<int> vertices;
                int leaf;
                std::vector<int> cubes;
                std::vector<std::vector<int>> embeddings;
            };
            std::vector<Component> components;

            // each overlap under the root, by its partition of the pattern vertices (sorted
            // blocks, sorted, singletons included)
            std::map<std::vector<std::vector<int>>, int> overlaps_by_partition;
        };

        // Delete everything but the counting formula from the core, which shows that the
        // root implies the encoding, and conclude that the core is equienumerable with
        // the translator's counting formula (phase 5 of decomposition-counting.md).
        auto finish(const Structure &) -> void;
    };
}

#endif
