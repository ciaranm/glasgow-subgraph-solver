/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_SYMMETRIES_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_SYMMETRIES_HH 1

#include "formats/input_graph.hh"

#include <exception>
#include <list>
#include <string>

class GapFailedUs :
    public std::exception
{
    private:
        std::string _what;

    public:
        GapFailedUs(const std::string & message) noexcept;

        auto what() const noexcept -> const char *;
};

auto find_symmetries(const char * const argv0, const InputGraph & graph, std::list<std::pair<std::string, std::string> > & constraints, std::string & aut_size) -> void;

#endif
