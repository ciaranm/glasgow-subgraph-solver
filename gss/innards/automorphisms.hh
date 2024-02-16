#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_GSS_INNARDS_AUTOMORPHISMS_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_GSS_INNARDS_AUTOMORPHISMS_HH

#include <gss/formats/input_graph.hh>
#include <gss/loooong.hh>
#include<gss/innards/svo_bitset.hh>

#include <list>
#include <string>
#include <utility>

namespace gss::innards
{
    using OrderConstraints = std::list<std::pair<std::string, std::string>>;

    auto automorphisms_as_order_constraints(const InputGraph &, const bool with_generators) -> OrderConstraints;
    auto automorphisms_as_order_constraints(const InputGraph &, std::vector<int> base) -> OrderConstraints;
    auto automorphisms_as_order_constraints_with_generators(const InputGraph &, std::vector<int> base) -> OrderConstraints;
    auto dynamic_order_constraints(std::vector<gss::innards::SVOBitset> m, std::vector<int> base) -> std::vector<std::pair<unsigned int, unsigned int>>;
}

#endif
