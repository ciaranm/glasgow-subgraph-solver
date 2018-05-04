/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef CODE_GUARD_RESULT_HH
#define CODE_GUARD_RESULT_HH 1

#include <map>
#include <list>
#include <chrono>

struct Result
{
    /// The isomorphism, empty if none found.
    std::map<int, int> isomorphism;

    /// Total number of nodes processed.
    unsigned long long nodes = 0;

    /// Extra stats, to output
    std::list<std::string> extra_stats;

    /// Number of solutions, only if enumerating
    unsigned long long solution_count = 0;
};

#endif
