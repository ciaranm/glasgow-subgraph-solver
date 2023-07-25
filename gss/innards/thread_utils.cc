#include <gss/innards/thread_utils.hh>

#include <thread>

using std::thread;

auto gss::innards::how_many_threads(unsigned n) -> unsigned
{
    if (0 == n)
        n = thread::hardware_concurrency();
    if (0 == n)
        n = 1;
    return n;
}
