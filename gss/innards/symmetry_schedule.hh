#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_SYMMETRY_SCHEDULE_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_SYMMETRY_SCHEDULE_HH 1

#include <chrono>

namespace gss::innards
{
    class SymmetrySchedule
    {
        public: 
            virtual ~SymmetrySchedule() = default;

            virtual auto did_a_backtrack() -> void = 0;
            virtual auto should_engage_symmetry() -> bool = 0;
    };

    /**
     * Do symmetry breaking from the beginning, if we're going to
     */
    class NoSymmetrySchedule final : public SymmetrySchedule 
    {
        public:
            NoSymmetrySchedule();

            virtual auto did_a_backtrack() -> void override;
            virtual auto should_engage_symmetry() -> bool override;
    };

    /**
     * Engage symmetry breaking after n microseconds of solve time
     */
    class TimedSymmetrySchedule final : public SymmetrySchedule
    {
        private:
            std::chrono::microseconds _duration;
            std::chrono::steady_clock::time_point _start_here;
        public:
            static constexpr std::chrono::microseconds default_duration{100000};

            TimedSymmetrySchedule(std::chrono::microseconds duration);

            virtual auto did_a_backtrack() -> void override;
            virtual auto should_engage_symmetry() -> bool override;
    };

    /**
     * Engage symmetry breaking after n backtracks
     */
    class BacktrackSymmetrySchedule final : public SymmetrySchedule
    {
        private:
            long long _number_of_backtracks = 0, _minimum_backtracks;
        public:
            static constexpr long long default_minimum_backtracks = 1;

            BacktrackSymmetrySchedule(long long minimum_backtracks);

            virtual auto did_a_backtrack() -> void override;
            virtual auto should_engage_symmetry() -> bool override;

    };
}

#endif