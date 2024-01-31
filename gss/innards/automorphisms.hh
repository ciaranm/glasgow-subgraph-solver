#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_GSS_INNARDS_AUTOMORPHISMS_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_GSS_INNARDS_AUTOMORPHISMS_HH

#include <gss/formats/input_graph.hh>
#include <gss/loooong.hh>

#include <list>
#include <string>
#include <utility>

//TODO include dejavu in here somehow?
// #include <dejavu.h>

namespace gss::innards
{
    using OrderConstraints = std::list<std::pair<std::string, std::string>>;

    auto automorphisms_as_order_constraints(const InputGraph &) -> OrderConstraints;
    auto automorphisms_as_order_constraints(const InputGraph &, std::vector<int> base) -> OrderConstraints;
    // auto dynamic_order_constraints(dejavu::static_graph g, std::vector<int> base, int size) -> std::vector<std::pair<unsigned int, unsigned int>>;
}

#endif
