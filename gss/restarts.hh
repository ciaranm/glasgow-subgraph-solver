#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_RESTARTS_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_RESTARTS_HH 1

#include <atomic>
#include <chrono>
#include <list>
#include <memory>

namespace gss
{
    class RestartsSchedule
    {
    public:
        virtual ~RestartsSchedule() = default;

        virtual auto did_a_backtrack() -> void = 0;
        virtual auto did_a_restart() -> void = 0;
        virtual auto should_restart() -> bool = 0;
        virtual auto might_restart() -> bool = 0;
        virtual auto clone() -> RestartsSchedule * = 0;
    };

    class NoRestartsSchedule final : public RestartsSchedule
    {
    public:
        virtual auto did_a_backtrack() -> void override;
        virtual auto did_a_restart() -> void override;
        virtual auto should_restart() -> bool override;
        virtual auto might_restart() -> bool override;
        virtual auto clone() -> NoRestartsSchedule * override;
    };

    class LubyRestartsSchedule final : public RestartsSchedule
    {
    private:
        long long _backtracks_remaining;
        std::list<long long> _sequence;
        std::list<long long>::const_iterator _current_sequence;

    public:
        static constexpr unsigned long long default_multiplier = 666; // chosen by divine inspiration

        explicit LubyRestartsSchedule(long long multiplier);
        LubyRestartsSchedule(const LubyRestartsSchedule &);

        virtual auto did_a_backtrack() -> void override;
        virtual auto did_a_restart() -> void override;
        virtual auto should_restart() -> bool override;
        virtual auto might_restart() -> bool override;
        virtual auto clone() -> LubyRestartsSchedule * override;
    };

    class GeometricRestartsSchedule final : public RestartsSchedule
    {
    private:
        long long _number_of_backtracks = 0;
        double _current_value, _multiplier;

    public:
        static constexpr double default_initial_value = 5400;
        static constexpr double default_multiplier = 1.0;

        GeometricRestartsSchedule(double initial_value, double multiplier);

        virtual auto did_a_backtrack() -> void override;
        virtual auto did_a_restart() -> void override;
        virtual auto should_restart() -> bool override;
        virtual auto might_restart() -> bool override;
        virtual auto clone() -> GeometricRestartsSchedule * override;
    };

    class SyncedRestartSchedule final : public RestartsSchedule
    {
    private:
        std::atomic<bool> & _synchroniser;

    public:
        explicit SyncedRestartSchedule(std::atomic<bool> &);

        virtual auto did_a_backtrack() -> void override;
        virtual auto did_a_restart() -> void override;
        virtual auto should_restart() -> bool override;
        virtual auto might_restart() -> bool override;
        virtual auto clone() -> SyncedRestartSchedule * override;
    };

    class TimedRestartsSchedule final : public RestartsSchedule
    {
    private:
        long long _number_of_backtracks = 0, _minimum_backtracks;
        std::chrono::milliseconds _duration;
        std::chrono::steady_clock::time_point _next_restart_point;

    public:
        static constexpr std::chrono::milliseconds default_duration{100};
        static constexpr unsigned long long default_minimum_backtracks = 100;

        TimedRestartsSchedule(std::chrono::milliseconds duration, unsigned long long minimum_backtracks);

        virtual auto did_a_backtrack() -> void override;
        virtual auto did_a_restart() -> void override;
        virtual auto should_restart() -> bool override;
        virtual auto might_restart() -> bool override;
        virtual auto clone() -> TimedRestartsSchedule * override;
    };
}

#endif
