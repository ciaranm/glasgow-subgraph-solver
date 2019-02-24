/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "restarts.hh"
#include <algorithm>
#include <cmath>

using std::advance;
using std::distance;
using std::round;

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

