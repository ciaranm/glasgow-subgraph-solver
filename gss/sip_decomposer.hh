#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_SIP_DECOMPOSER_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_SIP_DECOMPOSER_HH 1

#include <gss/homomorphism.hh>

#include <iosfwd>
#include <optional>

namespace gss
{
    struct DecompositionOptions
    {
        /// If not null, also write a counting statement for the count here (see
        /// dev_docs/counting-statements.md), unless the count is zero because nothing can
        /// embed, in which case nothing is written. Writing a statement is only supported
        /// when counting injective, non-induced embeddings of an undirected pattern.
        std::ostream * statement = nullptr;

        /// How many of the target's highest-degree vertices to treat as hubs, onto which
        /// pattern vertices are guessed rather than searched for (c in Section 3.6 of
        /// Jessica Ryan's thesis). If not given, this is chosen from the target's degrees.
        /// Any choice gives the same count. Only used when counting; with proof logging,
        /// there are no hubs unless this says otherwise, which is not yet supported.
        std::optional<int> hubs;

        /// Whether to start with the plain solver's checks before search, on the whole
        /// instance, which can rule out every embedding where each component would fit on
        /// its own. Not yet with proof logging. Only worth turning off to exercise the
        /// decomposition itself, as the tests do.
        bool check_whole_instance = true;
    };

    auto solve_sip_by_decomposition(
        const InputGraph & pattern,
        const InputGraph & target,
        const HomomorphismParams & params,
        const DecompositionOptions & options = {}) -> HomomorphismResult;
}

#endif
