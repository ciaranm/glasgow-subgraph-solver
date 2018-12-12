/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_RESTARTS_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_RESTARTS_HH 1

#include <list>

class RestartsSchedule
{
    public:
        virtual ~RestartsSchedule() = default;

        virtual auto did_a_backtrack() -> void = 0;
        virtual auto did_a_restart() -> void = 0;
        virtual auto should_restart() -> bool = 0;
};

class NoRestartsSchedule final : public RestartsSchedule
{
    public:
        virtual auto did_a_backtrack() -> void override;
        virtual auto did_a_restart() -> void override;
        virtual auto should_restart() -> bool override;
};

class LubyRestartsSchedule final : public RestartsSchedule
{
    private:
        long long _backtracks_remaining;
        std::list<long long> _sequence;
        std::list<long long>::const_iterator _current_sequence;

    public:
        static constexpr unsigned long long default_multiplier = 660; // chosen by SMAC

        explicit LubyRestartsSchedule(long long multiplier);

        virtual auto did_a_backtrack() -> void override;
        virtual auto did_a_restart() -> void override;
        virtual auto should_restart() -> bool override;
};

class GeometricRestartsSchedule final : public RestartsSchedule
{
    private:
        long long _number_of_backtracks = 0;
        double _current_value, _multiplier;

    public:
        static constexpr double default_initial_value = 660;
        static constexpr double default_multiplier = 1.0;

        GeometricRestartsSchedule(double initial_value, double multiplier);

        virtual auto did_a_backtrack() -> void override;
        virtual auto did_a_restart() -> void override;
        virtual auto should_restart() -> bool override;
};

#endif
