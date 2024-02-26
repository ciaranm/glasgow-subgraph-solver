#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_GSS_INNARDS_AUTOMORPHISMS_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_GSS_INNARDS_AUTOMORPHISMS_HH

#include <gss/formats/input_graph.hh>
#include <gss/loooong.hh>
#include <gss/innards/svo_bitset.hh>

#include <list>
#include <string>
#include <utility>

#include "build/_deps/dejavu-src/dejavu.h"

namespace gss::innards
{
    using OrderConstraints = std::list<std::pair<std::string, std::string>>;

    auto automorphisms_as_order_constraints(const InputGraph &, const bool with_generators) -> OrderConstraints;
    auto automorphisms_as_order_constraints(const InputGraph &, std::vector<int> base) -> OrderConstraints;
    auto automorphisms_as_order_constraints_with_generators(const InputGraph &, std::vector<int> base) -> OrderConstraints;
    auto initialise_dynamic_structure(dejavu::groups::random_schreier &, std::vector<innards::SVOBitset> m) -> void;
    auto dynamic_order_constraints(int sz, std::vector<int> base, dejavu::groups::random_schreier &, std::vector<std::pair<unsigned int, unsigned int>> &) -> void;
}

#endif
