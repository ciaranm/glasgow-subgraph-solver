#include <gss/timeout.hh>

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>

#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>

using namespace gss;

using std::atomic;
using std::condition_variable;
using std::cv_status;
using std::make_unique;
using std::mutex;
using std::thread;
using std::unique_lock;

using std::chrono::operator""s;
using std::chrono::seconds;
using std::chrono::system_clock;

struct Timeout::Detail
{
    bool aborted = false;
    thread timeout_thread;
    mutex timeout_mutex;
    condition_variable timeout_cv;
    atomic<bool> abort;
};

Timeout::Timeout(const seconds limit) :
    _detail(make_unique<Detail>())
{
    _detail->abort.store(false);
    if (0s != limit) {
        _detail->timeout_thread = thread([limit, &detail = this->_detail] {
            auto abort_time = system_clock::now() + limit;
            {
                /* Sleep until either we've reached the time limit,
                 * or we've finished all the work. */
                unique_lock<mutex> guard(detail->timeout_mutex);
                while (! detail->abort.load()) {
                    if (cv_status::timeout == detail->timeout_cv.wait_until(guard, abort_time)) {
                        /* We've woken up, and it's due to a timeout. */
                        detail->aborted = true;
                        break;
                    }
                }
            }
            detail->abort.store(true);
        });
    }
}

Timeout::~Timeout()
{
    stop();
}

auto Timeout::should_abort() const -> bool
{
    return _detail->abort.load();
}

auto Timeout::aborted() const -> bool
{
    return _detail->aborted;
}

auto Timeout::trigger_early_abort() -> void
{
    return _detail->abort.store(true);
}

auto Timeout::stop() -> void
{
    /* Clean up the timeout thread */
    if (_detail->timeout_thread.joinable()) {
        {
            unique_lock<mutex> guard(_detail->timeout_mutex);
            _detail->abort.store(true);
            _detail->timeout_cv.notify_all();
        }
        _detail->timeout_thread.join();
    }
}
