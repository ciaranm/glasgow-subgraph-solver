#include <gss/restarts.hh>

#include <algorithm>
#include <cmath>

using namespace gss;

using std::advance;
using std::distance;
using std::round;

using std::chrono::milliseconds;
using std::chrono::steady_clock;

auto NoRestartsSchedule::did_a_backtrack() -> void
{
}

auto NoRestartsSchedule::did_a_restart() -> void
{
}

auto NoRestartsSchedule::should_restart() -> bool
{
    return false;
}

auto NoRestartsSchedule::might_restart() -> bool
{
    return false;
}

auto NoRestartsSchedule::clone() -> NoRestartsSchedule *
{
    return new NoRestartsSchedule(*this);
}

LubyRestartsSchedule::LubyRestartsSchedule(long long m) :
    _backtracks_remaining(m)
{
    _current_sequence = _sequence.insert(_sequence.end(), m);
}

LubyRestartsSchedule::LubyRestartsSchedule(const LubyRestartsSchedule & other) :
    _backtracks_remaining(other._backtracks_remaining),
    _sequence(other._sequence),
    _current_sequence(_sequence.begin())
{
    advance(_current_sequence, distance(other._sequence.begin(), other._current_sequence));
}

auto LubyRestartsSchedule::did_a_backtrack() -> void
{
    --_backtracks_remaining;
}

auto LubyRestartsSchedule::did_a_restart() -> void
{
    if (next(_current_sequence) == _sequence.end()) {
        _sequence.insert(_sequence.end(), _sequence.begin(), _sequence.end());
        _sequence.push_back(_sequence.back() * 2);
    }
    ++_current_sequence;
    _backtracks_remaining = *_current_sequence;
}

auto LubyRestartsSchedule::should_restart() -> bool
{
    return _backtracks_remaining <= 0;
}

auto LubyRestartsSchedule::clone() -> LubyRestartsSchedule *
{
    return new LubyRestartsSchedule(*this);
}

auto LubyRestartsSchedule::might_restart() -> bool
{
    return true;
}

GeometricRestartsSchedule::GeometricRestartsSchedule(double v, double m) :
    _current_value(v),
    _multiplier(m)
{
}

auto GeometricRestartsSchedule::did_a_backtrack() -> void
{
    ++_number_of_backtracks;
}

auto GeometricRestartsSchedule::did_a_restart() -> void
{
    _number_of_backtracks = 0;
    _current_value *= _multiplier;
}

auto GeometricRestartsSchedule::should_restart() -> bool
{
    return _number_of_backtracks >= round(_current_value);
}

auto GeometricRestartsSchedule::clone() -> GeometricRestartsSchedule *
{
    return new GeometricRestartsSchedule(*this);
}

auto GeometricRestartsSchedule::might_restart() -> bool
{
    return true;
}

SyncedRestartSchedule::SyncedRestartSchedule(std::atomic<bool> & a) :
    _synchroniser(a)
{
}

auto SyncedRestartSchedule::did_a_backtrack() -> void
{
}

auto SyncedRestartSchedule::did_a_restart() -> void
{
}

auto SyncedRestartSchedule::should_restart() -> bool
{
    return _synchroniser;
}

auto SyncedRestartSchedule::might_restart() -> bool
{
    return true;
}

auto SyncedRestartSchedule::clone() -> SyncedRestartSchedule *
{
    return new SyncedRestartSchedule(*this);
}

TimedRestartsSchedule::TimedRestartsSchedule(milliseconds d, unsigned long long m) :
    _number_of_backtracks(0),
    _minimum_backtracks(m),
    _duration(d),
    _next_restart_point(steady_clock::now() + d)
{
}

auto TimedRestartsSchedule::did_a_backtrack() -> void
{
    ++_number_of_backtracks;
}

auto TimedRestartsSchedule::did_a_restart() -> void
{
    _next_restart_point = steady_clock::now() + _duration;
    _number_of_backtracks = 0;
}

auto TimedRestartsSchedule::should_restart() -> bool
{
    return (_number_of_backtracks >= _minimum_backtracks) && (steady_clock::now() >= _next_restart_point);
}

auto TimedRestartsSchedule::clone() -> TimedRestartsSchedule *
{
    return new TimedRestartsSchedule(*this);
}

auto TimedRestartsSchedule::might_restart() -> bool
{
    return true;
}
