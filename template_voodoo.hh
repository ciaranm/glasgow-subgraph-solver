/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_TEMPLATE_VOODOO_HH
#define GLASGOW_SUBGRAPH_SOLVER_TEMPLATE_VOODOO_HH 1

#include "fixed_bit_set.hh"

#include <array>
#include <type_traits>
#include <utility>
#include <vector>

#include <boost/dynamic_bitset.hpp>

template <template <typename, typename> class Algorithm_, typename Result_, typename Graph_, unsigned size_, unsigned... other_sizes_, typename... Params_>
auto select_graph_size(const std::integer_sequence<unsigned, size_, other_sizes_...> &, const Graph_ & graph, Params_ && ... params) -> Result_
{
    if (graph.size() < int(size_ * bits_per_word)) {
        Algorithm_<FixedBitSet<size_>, std::array<int, size_ * bits_per_word + 1> > algorithm{ graph, std::forward<Params_>(params)... };
        return algorithm.run();
    }
    else {
        if constexpr (0 == sizeof...(other_sizes_)) {
            Algorithm_<boost::dynamic_bitset<>, std::vector<int> > algorithm{ graph, std::forward<Params_>(params)... };
            return algorithm.run();
        }
        else
            return select_graph_size<Algorithm_, Result_, Graph_>(std::integer_sequence<unsigned, other_sizes_...>{}, graph, std::forward<Params_>(params)...);
    }
}

using AllGraphSizes = std::integer_sequence<unsigned, 1, 2, 3, 4, 5, 6, 7, 8, 16, 20, 24, 28, 32, 64, 128, 256, 512, 1024>;

#endif
