#include <gss/restarts.hh>

#include <catch2/catch_test_macros.hpp>

#include <atomic>
#include <memory>
#include <vector>

using namespace gss;

using std::atomic;
using std::unique_ptr;
using std::vector;

using std::chrono::hours;
using std::chrono::milliseconds;

namespace
{
    // Run a schedule forward one restart's worth of backtracks, returning how many
    // backtracks it took to trigger the restart.
    auto backtracks_until_restart(RestartsSchedule & schedule) -> long long
    {
        long long n = 0;
        while (! schedule.should_restart()) {
            schedule.did_a_backtrack();
            ++n;
        }
        schedule.did_a_restart();
        return n;
    }

    auto first_restart_lengths(RestartsSchedule & schedule, unsigned how_many) -> vector<long long>
    {
        vector<long long> lengths;
        for (unsigned i = 0; i < how_many; ++i)
            lengths.push_back(backtracks_until_restart(schedule));
        return lengths;
    }
}

TEST_CASE("NoRestartsSchedule never restarts")
{
    NoRestartsSchedule schedule;
    CHECK_FALSE(schedule.might_restart());
    CHECK_FALSE(schedule.should_restart());

    for (int i = 0; i < 1000; ++i) {
        schedule.did_a_backtrack();
        CHECK_FALSE(schedule.should_restart());
    }
}

TEST_CASE("LubyRestartsSchedule follows the Luby sequence")
{
    // The (unscaled) Luby sequence.
    const vector<long long> luby{1, 1, 2, 1, 1, 2, 4, 1, 1, 2, 1, 1, 2, 4, 8};

    SECTION("multiplier 1")
    {
        LubyRestartsSchedule schedule{1};
        CHECK(schedule.might_restart());
        CHECK(first_restart_lengths(schedule, luby.size()) == luby);
    }

    SECTION("scaled by a multiplier")
    {
        const long long multiplier = 7;
        LubyRestartsSchedule schedule{multiplier};
        vector<long long> expected;
        for (auto v : luby)
            expected.push_back(v * multiplier);
        CHECK(first_restart_lengths(schedule, luby.size()) == expected);
    }
}

TEST_CASE("GeometricRestartsSchedule grows by its multiplier")
{
    SECTION("doubling")
    {
        GeometricRestartsSchedule schedule{10.0, 2.0};
        CHECK(schedule.might_restart());
        CHECK(first_restart_lengths(schedule, 4) == vector<long long>{10, 20, 40, 80});
    }

    SECTION("constant when the multiplier is 1")
    {
        GeometricRestartsSchedule schedule{15.0, 1.0};
        CHECK(first_restart_lengths(schedule, 4) == vector<long long>{15, 15, 15, 15});
    }
}

TEST_CASE("SyncedRestartSchedule mirrors its synchroniser")
{
    atomic<bool> trigger{false};
    SyncedRestartSchedule schedule{trigger};

    CHECK(schedule.might_restart());
    CHECK_FALSE(schedule.should_restart());

    trigger = true;
    CHECK(schedule.should_restart());

    trigger = false;
    CHECK_FALSE(schedule.should_restart());
}

TEST_CASE("TimedRestartsSchedule respects the minimum backtracks / time gate")
{
    // With a duration of an hour, the time condition cannot be met during the test,
    // so the schedule must never ask to restart however many backtracks happen.
    TimedRestartsSchedule schedule{hours{1}, 10};
    CHECK(schedule.might_restart());

    for (int i = 0; i < 100; ++i) {
        schedule.did_a_backtrack();
        CHECK_FALSE(schedule.should_restart());
    }
}

TEST_CASE("RestartsSchedule clone produces an independent equivalent schedule")
{
    SECTION("Luby")
    {
        LubyRestartsSchedule original{1};
        // Advance the original past its first two restarts.
        backtracks_until_restart(original);
        backtracks_until_restart(original);

        unique_ptr<RestartsSchedule> copy{original.clone()};

        // The clone resumes from the same point in the sequence as the original.
        CHECK(backtracks_until_restart(*copy) == backtracks_until_restart(original));
        CHECK(backtracks_until_restart(*copy) == backtracks_until_restart(original));
    }

    SECTION("None")
    {
        NoRestartsSchedule original;
        unique_ptr<RestartsSchedule> copy{original.clone()};
        CHECK_FALSE(copy->should_restart());
    }
}
