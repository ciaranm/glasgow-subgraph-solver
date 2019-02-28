/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "thread_utils.hh"

#include <thread>

using std::thread;

auto how_many_threads(unsigned n) -> unsigned
{
    if (0 == n)
        n = thread::hardware_concurrency();
    if (0 == n)
        n = 1;
    return n;
}

