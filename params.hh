/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_PARAMS_HH
#define GLASGOW_SUBGRAPH_SOLVER_PARAMS_HH 1

#include <chrono>
#include <atomic>
#include <cmath>

struct Params
{
    /// If this is set to true, we should abort due to a time limit.
    std::atomic<bool> * abort;

    /// The start time of the algorithm.
    std::chrono::time_point<std::chrono::steady_clock> start_time;

    /// Induced?
    bool induced = false;

    /// Enumerate?
    bool enumerate = false;

    /// Presolve?
    bool presolve = false;

    /// Default chosen by SMAC
    static constexpr unsigned long long dodgy_default_magic_luby_multiplier = 660;

    /// Multiplier for Luby sequence
    unsigned long long luby_multiplier = dodgy_default_magic_luby_multiplier;
};

#endif
