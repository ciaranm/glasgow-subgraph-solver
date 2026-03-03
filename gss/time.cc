#include <gss/time.hh>

using namespace gss;

using std::string;
using std::to_string;
using std::chrono::microseconds;

auto gss::microseconds_to_string(const microseconds & v) -> string
{
    return to_string(v.count() / 1.0e6) + "s";
}
