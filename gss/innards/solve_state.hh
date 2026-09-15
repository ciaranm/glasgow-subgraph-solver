#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_SOLVE_STATE_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_SOLVE_STATE_HH 1

#include <gss/innards/homomorphism_domain.hh>
#include <gss/innards/homomorphism_model.hh>
#include <gss/innards/homomorphism_searcher.hh>

#include <memory>
#include <vector>

namespace gss::innards
{
    // The working state carried through the solve pipeline (see
    // dev_docs/preprocessor-refactor.md). It gathers the three things that were
    // previously locals scattered across the solver: the constraint model, the
    // (root) domains, and the nogood store. The model is "frozen" while a search
    // step runs and may be mutated between steps; the domains and nogoods persist
    // across steps. Folding them here is the prerequisite for staged solving
    // (Phase 6), where builder / filter steps grow the model and bounded-search
    // steps accumulate nogoods between rounds.
    //
    // The model is built lazily by the main solve step (after the cheap shortcut
    // steps have had their chance to conclude), so it starts null. The threaded
    // search keeps its own per-thread domain copies and per-thread nogood stores
    // and reads only the shared model from here; the carried domains / watches are
    // used by the sequential path (the one staged solving will build on).
    struct SolveState
    {
        std::unique_ptr<HomomorphismModel> model;
        std::vector<HomomorphismDomain> domains;
        Watches<HomomorphismAssignment, HomomorphismAssignmentWatchTable> watches;
    };
}

#endif
