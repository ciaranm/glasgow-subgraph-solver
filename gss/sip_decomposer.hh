#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_SIP_DECOMPOSER_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_SIP_DECOMPOSER_HH 1

#include <gss/counting_statement.hh>
#include <gss/homomorphism.hh>

namespace gss
{
    /// If statement is given, also writes into it a counting-language statement of the
    /// argument, whose count is the solution count (see
    /// dev_docs/decomposition-counting.md). That needs injective, non-induced,
    /// undirected counting, and is much slower than counting alone.
    auto solve_sip_by_decomposition(
        const InputGraph & pattern,
        const InputGraph & target,
        const HomomorphismParams & params,
        CountingStatement * statement = nullptr) -> HomomorphismResult;
}

#endif
