/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_PARAMS_HH
#define GLASGOW_SUBGRAPH_SOLVER_PARAMS_HH 1

#include <atomic>
#include <chrono>
#include <cmath>
#include <exception>
#include <limits>
#include <memory>
#include <string>

#include "restarts.hh"
#include "value_ordering.hh"

struct Params
{
    /// If this is set to true, we should abort due to a time limit.
    std::atomic<bool> * abort;

    /// The start time of the algorithm.
    std::chrono::time_point<std::chrono::steady_clock> start_time;

    /// Induced?
    bool induced = false;

    /// Noninjective?
    bool noninjective = false;

    /// Enumerate?
    bool enumerate = false;

    /// Which value-ordering heuristic?
    ValueOrdering value_ordering_heuristic = ValueOrdering::Biased;

    /// Restarts schedule
    std::unique_ptr<RestartsSchedule> restarts_schedule;

    /// Largest size of nogood to store (0 disables nogoods)
    unsigned nogood_size_limit = std::numeric_limits<unsigned>::max();
};

class UnsupportedConfiguration :
    public std::exception
{
    private:
        std::string _what;

    public:
        UnsupportedConfiguration(const std::string & message) throw ();

        auto what() const throw () -> const char *;
};

#endif
