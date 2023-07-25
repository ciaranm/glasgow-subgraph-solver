#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_LACKEY_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_LACKEY_HH 1

#include <gss/formats/input_graph.hh>
#include <gss/vertex_to_vertex_mapping.hh>

#include <exception>
#include <functional>
#include <memory>
#include <string>

namespace gss::innards
{
    class DisobedientLackeyError : public std::exception
    {
    private:
        std::string _what;

    public:
        explicit DisobedientLackeyError(const std::string & message) noexcept;

        auto what() const throw() -> const char *;
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
        Lackey & operator=(const Lackey &) = delete;

        using DeletionFunction = std::function<auto(int, int)->bool>;
        using RestrictRangeFunction = std::function<auto(int, int)->void>;

        auto check_solution(
            const VertexToVertexMapping &,
            bool partial,
            bool all_solutions,
            const DeletionFunction & deletions) -> bool;

        auto reduce_initial_bounds(
            const RestrictRangeFunction & restrict_range) -> bool;

        auto number_of_checks() const -> long;
        auto number_of_propagations() const -> long;
        auto number_of_deletions() const -> long;
        auto number_of_calls() const -> long;
    };
}

#endif
