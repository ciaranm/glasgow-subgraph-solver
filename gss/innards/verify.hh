#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_VERIFY_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_VERIFY_HH 1

#include <gss/formats/input_graph.hh>

#include <exception>
#include <map>
#include <string>
#include <utility>

namespace gss::innards
{
    auto verify_homomorphism(
        const InputGraph & pattern,
        const InputGraph & target,
        bool injective,
        bool locally_injective,
        bool induced,
        const std::map<int, int> & mapping) -> void;

    /**
     * The cost of a mapping, computed directly from the graphs rather than from anything
     * the solver built: the costs of the target vertices used, plus the cost of the
     * target edge each pattern edge lands on. A pattern edge is counted once per arc if
     * either graph is directed (so an undirected pattern edge in a directed target uses
     * both arcs), and once otherwise, which is what reification does.
     *
     * \throw BuggySolution if a pattern edge lands on no target edge with its label.
     */
    auto cost_of_mapping(
        const InputGraph & pattern,
        const InputGraph & target,
        const std::map<int, int> & mapping) -> long long;

    class BuggySolution : public std::exception
    {
    private:
        std::string _what;

    public:
        BuggySolution(const std::string & message) noexcept;

        auto what() const noexcept -> const char * override;
    };
}

#endif
