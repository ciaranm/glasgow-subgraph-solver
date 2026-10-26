#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_SIP_DECOMPOSER_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_SIP_DECOMPOSER_HH 1

#include <gss/homomorphism.hh>

#include <iosfwd>

namespace gss
{
    /// If statement is not null, also writes a counting statement for the count to it
    /// (see dev_docs/counting-statements.md), unless some component of the pattern has
    /// no embeddings at all, in which case nothing is written. Writing a statement is
    /// only supported when counting injective, non-induced embeddings of an undirected
    /// pattern.
    auto solve_sip_by_decomposition(
        const InputGraph & pattern,
        const InputGraph & target,
        const HomomorphismParams & params,
        std::ostream * statement = nullptr) -> HomomorphismResult;
}

#endif
