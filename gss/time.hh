#ifndef GLASGOW_SUBGRAPH_SOLVER_TIME_HH
#define GLASGOW_SUBGRAPH_SOLVER_TIME_HH

#include <chrono>
#include <string>

namespace gss
{
    [[nodiscard]] auto microseconds_to_string(const std::chrono::microseconds & v) -> std::string;
}

#endif
