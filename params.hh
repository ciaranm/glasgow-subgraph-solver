/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef CODE_GUARD_PARAMS_HH
#define CODE_GUARD_PARAMS_HH 1

#include <chrono>
#include <atomic>
#include <cmath>

struct Params
{
    /// If this is set to true, we should abort due to a time limit.
    std::atomic<bool> * abort;

    /// The start time of the algorithm.
    std::chrono::time_point<std::chrono::steady_clock> start_time;

    /// Use dds?
    bool dds = false;

    /// Use restarts?
    bool restarts = false;

    /// Use no heuristic?
    bool input_order = false;

    /// Use shuffles? (For science purposes, not real use)
    bool shuffle = false;

    /// Use softmax shuffles?
    bool softmax_shuffle = false;

    /// Use antiheuristic? (For science purposes)
    bool antiheuristic = false;

    /// Don't use nogoods?
    bool goods = false;

    /// Induced?
    bool induced = false;

    /// Enumerate?
    bool enumerate = false;

    /// Default chosen by SMAC
    static constexpr unsigned long long dodgy_default_magic_luby_multiplier = 660;

    /// Multiplier for Luby sequence
    unsigned long long luby_multiplier = dodgy_default_magic_luby_multiplier;

    /// Multiplier for geometric sequence (set to 0 for luby)
    unsigned long long geometric_multiplier = 0.0;

    /// Initial geometric sequence value
    unsigned long long geometric_start = 10;

    /// Specify a random seed.
    unsigned seed = 0;
};

#endif
