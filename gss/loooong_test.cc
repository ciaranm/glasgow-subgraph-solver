#include <gss/loooong.hh>

#include <catch2/catch_test_macros.hpp>

#include <sstream>
#include <string>
#include <utility>

using namespace gss;

using std::string;
using std::stringstream;

namespace
{
    auto to_string(const loooong & x) -> string
    {
        stringstream s;
        s << x;
        return s.str();
    }
}

TEST_CASE("loooong construction, equality and printing")
{
    CHECK(loooong{42} == loooong{42L});
    CHECK(loooong{42} == loooong{42UL});
    CHECK(loooong{-5} != loooong{5});

    CHECK(to_string(loooong{0}) == "0");
    CHECK(to_string(loooong{12345}) == "12345");
    CHECK(to_string(loooong{-12345}) == "-12345");
}

TEST_CASE("loooong arithmetic operators")
{
    const loooong a{20}, b{7};
    CHECK(a + b == loooong{27});
    CHECK(a - b == loooong{13});
    CHECK(a * b == loooong{140});
    CHECK(a / b == loooong{2}); // truncating division
    CHECK(loooong{3} - loooong{10} == loooong{-7});
}

TEST_CASE("loooong compound assignment operators")
{
    loooong c{10};
    c += loooong{5};
    CHECK(c == loooong{15});
    c -= loooong{20};
    CHECK(c == loooong{-5});
    c *= loooong{4};
    CHECK(c == loooong{-20});
    c /= loooong{3}; // -20 / 3 truncates toward zero
    CHECK(c == loooong{-6});
    ++c;
    CHECK(c == loooong{-5});
}

TEST_CASE("loooong comparisons")
{
    CHECK(loooong{3} < loooong{5});
    CHECK(loooong{5} > loooong{3});
    CHECK(loooong{5} <= loooong{5});
    CHECK(loooong{5} >= loooong{5});
    CHECK(loooong{-1} < loooong{0});
}

TEST_CASE("loooong handles values beyond 64 bits")
{
    // 2^64 = 18446744073709551616, one more bit than unsigned long long holds.
    loooong power{1};
    for (int i = 0; i < 64; ++i)
        power *= loooong{2};
    CHECK(to_string(power) == "18446744073709551616");

    // 25! = 15511210043330985984000000.
    loooong factorial{1};
    for (int i = 1; i <= 25; ++i)
        factorial *= loooong{i};
    CHECK(to_string(factorial) == "15511210043330985984000000");
    CHECK(factorial > loooong{0});
}

TEST_CASE("loooong gcd")
{
    CHECK(gcd(loooong{12}, loooong{18}) == loooong{6});
    CHECK(gcd(loooong{17}, loooong{5}) == loooong{1});
    CHECK(gcd(loooong{0}, loooong{9}) == loooong{9});
}

TEST_CASE("loooong copies and moves are independent")
{
    loooong a{100};

    loooong copy{a};
    copy += loooong{1};
    CHECK(a == loooong{100}); // original untouched
    CHECK(copy == loooong{101});

    loooong assigned;
    assigned = a;
    assigned += loooong{5};
    CHECK(a == loooong{100});
    CHECK(assigned == loooong{105});

    loooong moved{std::move(copy)};
    CHECK(moved == loooong{101});
}
