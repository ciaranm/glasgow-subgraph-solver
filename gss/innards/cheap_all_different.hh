#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_CHEAP_ALL_DIFFERENT_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_CHEAP_ALL_DIFFERENT_HH 1

#include <gss/innards/homomorphism_domain.hh>
#include <gss/innards/homomorphism_model.hh>
#include <gss/innards/proof.hh>

#include <optional>
#include <vector>

namespace gss::innards
{
    /**
     * \param only_below if given, consider only the domains of pattern vertices numbered
     *     below this, leaving the others alone. For a reified instance, where every
     *     edge-vertex is kept distinct by its endpoints being distinct (see
     *     reification.hh), that is the original vertices.
     */
    auto cheap_all_different(unsigned target_size, std::vector<HomomorphismDomain> & domains, const std::shared_ptr<Proof> & proof,
        const HomomorphismModel * const, std::optional<unsigned> only_below = std::nullopt) -> bool;
}

#endif
