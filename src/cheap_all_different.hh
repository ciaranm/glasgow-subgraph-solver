/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_CHEAP_ALL_DIFFERENT_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_CHEAP_ALL_DIFFERENT_HH 1

#include "homomorphism_domain.hh"
#include "homomorphism_model.hh"
#include "proof.hh"

#include <vector>

auto cheap_all_different(unsigned target_size, std::vector<HomomorphismDomain> & domains, const std::shared_ptr<Proof> & proof,
        const HomomorphismModel * const) -> bool;

#endif
