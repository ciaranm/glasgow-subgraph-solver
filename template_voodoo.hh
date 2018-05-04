/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef CODE_GUARD_TEMPLATE_VOODOO_HH
#define CODE_GUARD_TEMPLATE_VOODOO_HH 1

#include <type_traits>

template <unsigned...>
struct GraphSizes;

struct NoMoreGraphSizes
{
};

template <unsigned n_, unsigned... n_rest_>
struct GraphSizes<n_, n_rest_...>
{
    enum { n = n_ };

    using Rest = GraphSizes<n_rest_...>;
};

template <unsigned n_>
struct GraphSizes<n_>
{
    enum { n = n_ };

    using Rest = NoMoreGraphSizes;
};

template <template <unsigned> class Algorithm_, typename Result_, typename Graph_, unsigned... sizes_, typename... Params_>
auto select_graph_size(const GraphSizes<sizes_...> &, const Graph_ & graph, Params_ && ... params) -> Result_
{
    if (graph.size() < GraphSizes<sizes_...>::n * bits_per_word) {
        Algorithm_<GraphSizes<sizes_...>::n> algorithm{
            graph, std::forward<Params_>(params)... };
        return algorithm.run();
    }
    else
        return select_graph_size<Algorithm_, Result_, Graph_>(typename GraphSizes<sizes_...>::Rest(), graph, std::forward<Params_>(params)...);
}

template <template <unsigned> class Algorithm_, typename Result_, typename Graph_, typename... Params_>
auto select_graph_size(const NoMoreGraphSizes &, const Graph_ &, Params_ && ...) -> Result_
{
    throw GraphTooBig();
}

using AllGraphSizes = GraphSizes<1, 2, 3, 4, 5, 6, 7, 8, 16, 20, 24, 28, 32, 64, 128, 256, 512, 1024>;

#endif
