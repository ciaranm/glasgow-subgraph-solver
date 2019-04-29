/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_CONFIGURATION_HH
#define GLASGOW_SUBGRAPH_SOLVER_CONFIGURATION_HH 1

#include <exception>
#include <string>


class UnsupportedConfiguration :
    public std::exception
{
    private:
        std::string _what;

    public:
        UnsupportedConfiguration(const std::string & message) noexcept;

        auto what() const noexcept -> const char * override;
};

#endif
