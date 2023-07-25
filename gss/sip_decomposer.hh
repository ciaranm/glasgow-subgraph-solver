#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_SIP_DECOMPOSER_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_SIP_DECOMPOSER_HH 1

#include <gss/homomorphism.hh>

namespace gss
{
    auto solve_sip_by_decomposition(
        const InputGraph & pattern,
        const InputGraph & target,
        const HomomorphismParams & params) -> HomomorphismResult;
}

#endif
