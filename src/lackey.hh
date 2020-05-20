/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_LACKEY_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_LACKEY_HH 1

#include "formats/input_graph.hh"
#include "vertex_to_vertex_mapping.hh"

#include <exception>
#include <memory>
#include <string>

class DisobedientLackeyError :
    public std::exception
{
    private:
        std::string _what;

    public:
        explicit DisobedientLackeyError(const std::string & message) noexcept;

        auto what() const throw () -> const char *;
};

class Lackey
{
    private:
        struct Imp;
        std::unique_ptr<Imp> _imp;

    public:
        Lackey(
                const std::string & send_to_name,
                const std::string & read_from_name,
                const InputGraph & pattern,
                const InputGraph & target);
        ~Lackey();

        Lackey(const Lackey &) = delete;
        Lackey & operator= (const Lackey &) = delete;

        auto check_solution(const VertexToVertexMapping &, bool partial, bool all_solutions) -> bool;
};

#endif
