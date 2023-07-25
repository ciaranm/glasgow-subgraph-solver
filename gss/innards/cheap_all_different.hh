#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_CHEAP_ALL_DIFFERENT_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_CHEAP_ALL_DIFFERENT_HH 1

#include <gss/innards/homomorphism_domain.hh>
#include <gss/innards/homomorphism_model.hh>
#include <gss/innards/proof.hh>

#include <vector>

namespace gss::innards
{
    auto cheap_all_different(unsigned target_size, std::vector<HomomorphismDomain> & domains, const std::shared_ptr<Proof> & proof,
        const HomomorphismModel * const) -> bool;
}

#endif
