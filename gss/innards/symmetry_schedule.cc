#include <gss/innards/symmetry_schedule.hh>

using namespace gss;
using namespace gss::innards;

using std::chrono::microseconds;
using std::chrono::steady_clock;

NoSymmetrySchedule::NoSymmetrySchedule()
{    
}

auto NoSymmetrySchedule::did_a_backtrack() -> void 
{
}

auto NoSymmetrySchedule::should_engage_symmetry() -> bool 
{
    return true;
}

TimedSymmetrySchedule::TimedSymmetrySchedule(microseconds d) :
    _duration(d),
    _start_here(steady_clock::now() + d)
{    
}

auto TimedSymmetrySchedule::did_a_backtrack() -> void 
{
}

auto TimedSymmetrySchedule::should_engage_symmetry() -> bool 
{
    return steady_clock::now() >= _start_here;
}

BacktrackSymmetrySchedule::BacktrackSymmetrySchedule(long long n) :
    _number_of_backtracks(0),
    _minimum_backtracks(n)
{
}

auto BacktrackSymmetrySchedule::did_a_backtrack() -> void
{
    ++_number_of_backtracks;
}

auto BacktrackSymmetrySchedule::should_engage_symmetry() -> bool
{
    return _number_of_backtracks >= _minimum_backtracks;
}
