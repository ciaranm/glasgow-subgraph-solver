/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_DOMAIN_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_HOMOMORPHISM_DOMAIN_HH 1

#include "svo_bitset.hh"

struct HomomorphismDomain
{
    unsigned v;
    unsigned count;
    bool fixed = false;
    SVOBitset values;

    explicit HomomorphismDomain(unsigned s) :
        values(s, 0)
    {
    }

    HomomorphismDomain(const HomomorphismDomain &) = default;
    HomomorphismDomain(HomomorphismDomain &&) = default;
};

#endif
