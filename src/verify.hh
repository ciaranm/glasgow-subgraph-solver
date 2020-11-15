/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_VERIFY_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_VERIFY_HH 1

#include "formats/input_graph.hh"

#include <string>
#include <exception>
#include <map>
#include <utility>

auto verify_homomorphism(
        const InputGraph & pattern,
        const InputGraph & target,
        bool injective,
        bool locally_injective,
        bool induced,
        const std::map<int, int> & mapping) -> void;

class BuggySolution :
    public std::exception
{
    private:
        std::string _what;

    public:
        BuggySolution(const std::string & message) noexcept;

        auto what() const noexcept -> const char * override;
};

#endif
