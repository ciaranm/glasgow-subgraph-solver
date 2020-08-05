/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_SIP_DECOMPOSER_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_SIP_DECOMPOSER_HH 1

#include "homomorphism.hh"

auto solve_sip_by_decomposition(
        const InputGraph & pattern,
        const InputGraph & target,
        const HomomorphismParams & params) -> HomomorphismResult;

#endif
