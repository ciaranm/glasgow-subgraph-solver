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

    auto automorphisms_as_order_constraints(const InputGraph &, const bool with_generators, const bool degree_sequence, std::vector<int> &orbit_sizes, std::vector<int> &base) -> OrderConstraints;
    auto automorphisms_as_order_constraints(const InputGraph &, std::vector<int> & base, std::vector<int> & orbit_sizes, const bool degree_sequence) -> OrderConstraints;
    auto automorphisms_as_order_constraints_with_generators(const InputGraph &, std::vector<int> base, std::vector<int> & orbit_sizes) -> OrderConstraints;
    auto initialise_dynamic_structure(dejavu::groups::random_schreier &, std::vector<innards::SVOBitset> m, const bool directed) -> bool;
    auto dynamic_order_constraints(int sz, std::vector<int> & base, std::vector<int> & orbit_sizes, dejavu::groups::random_schreier &, std::vector<std::pair<unsigned int, unsigned int>> &) -> void;
    auto generating_set(const InputGraph &i) -> std::vector<std::vector<unsigned int>>;
    auto generating_set(const InputGraph &i, std::vector<int> base) -> std::vector<std::vector<unsigned int>>;
    auto dynamic_generating_set(std::vector<int> & base, int sz, dejavu::groups::random_schreier & rschreier, std::vector<std::vector<int>> & generators) -> void;
    auto invert_automorphism(std::vector<unsigned int> aut) -> std::vector<unsigned int>;
    auto invert_list(std::vector<std::vector<unsigned int>> ls) -> std::vector<std::vector<unsigned int>>;

}

#endif
