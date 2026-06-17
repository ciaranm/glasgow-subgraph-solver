#include <gss/innards/svo_bitset.hh>

#include <catch2/catch_test_macros.hpp>

using namespace gss::innards;

namespace
{
    // SVOBitset stores up to 16 64-bit words inline ("short") and spills to the heap
    // ("long") beyond that. Every behavioural test is run at a small size (<= 1024
    // bits, inline) and a large size (> 1024 bits, heap) so both representations and
    // the boundary between them are exercised.
    constexpr unsigned small_bits = 128;  // 2 words, inline
    constexpr unsigned large_bits = 2048; // 32 words, heap

    auto check_basic_ops(unsigned size) -> void
    {
        SVOBitset b(size, 0);
        CHECK_FALSE(b.any());
        CHECK(b.count() == 0);
        CHECK(b.find_first() == SVOBitset::npos);

        b.set(5);
        b.set(63);
        b.set(64);
        b.set(100);
        CHECK(b.any());
        CHECK(b.count() == 4);
        CHECK(b.test(5));
        CHECK(b.test(64));
        CHECK_FALSE(b.test(6));
        CHECK(b.find_first() == 5);

        b.reset(5);
        CHECK_FALSE(b.test(5));
        CHECK(b.count() == 3);
        CHECK(b.find_first() == 63);

        b.reset();
        CHECK_FALSE(b.any());
        CHECK(b.count() == 0);
        CHECK(b.find_first() == SVOBitset::npos);
    }

    auto check_bitwise_ops(unsigned size) -> void
    {
        SVOBitset a(size, 0), c(size, 0);
        a.set(1);
        a.set(2);
        a.set(3);
        c.set(2);
        c.set(3);
        c.set(4);

        SVOBitset intersection = a;
        intersection &= c;
        CHECK(intersection.count() == 2);
        CHECK(intersection.test(2));
        CHECK(intersection.test(3));
        CHECK_FALSE(intersection.test(1));

        SVOBitset uni = a;
        uni |= c;
        CHECK(uni.count() == 4);
        CHECK(uni.test(1));
        CHECK(uni.test(4));

        SVOBitset difference = a;
        difference.intersect_with_complement(c);
        CHECK(difference.count() == 1);
        CHECK(difference.test(1));
        CHECK_FALSE(difference.test(2));
    }

    auto check_copy_independence(unsigned size) -> void
    {
        SVOBitset a(size, 0);
        a.set(10);
        a.set(20);

        SVOBitset copy_constructed = a;
        copy_constructed.set(30);
        CHECK(a.count() == 2); // original untouched
        CHECK_FALSE(a.test(30));
        CHECK(copy_constructed.count() == 3);

        SVOBitset copy_assigned(size, 0);
        copy_assigned = a;
        copy_assigned.reset(10);
        CHECK(a.test(10)); // original untouched
        CHECK_FALSE(copy_assigned.test(10));
        CHECK(copy_assigned.test(20));
    }
}

TEST_CASE("SVOBitset basic operations, inline")
{
    check_basic_ops(small_bits);
}

TEST_CASE("SVOBitset basic operations, heap")
{
    check_basic_ops(large_bits);
}

TEST_CASE("SVOBitset bitwise operations, inline")
{
    check_bitwise_ops(small_bits);
}

TEST_CASE("SVOBitset bitwise operations, heap")
{
    check_bitwise_ops(large_bits);
}

TEST_CASE("SVOBitset copies are independent, inline")
{
    check_copy_independence(small_bits);
}

TEST_CASE("SVOBitset copies are independent, heap")
{
    check_copy_independence(large_bits);
}

TEST_CASE("SVOBitset addresses bits in the heap region")
{
    SVOBitset b(large_bits, 0);
    b.set(1025);
    b.set(2000);
    CHECK(b.test(1025));
    CHECK(b.test(2000));
    CHECK(b.count() == 2);
    CHECK(b.find_first() == 1025);
}

TEST_CASE("SVOBitset copy-assignment across the inline/heap boundary")
{
    SVOBitset small(small_bits, 0);
    small.set(5);
    SVOBitset large(large_bits, 0);
    large.set(1500);

    SVOBitset target(small_bits, 0);
    target = large; // inline -> heap (allocates)
    CHECK(target.test(1500));
    CHECK(target.count() == 1);

    target = small; // heap -> inline (frees)
    CHECK(target.test(5));
    // target is now a 128-bit bitset again, so bit 1500 is out of range and must
    // not be probed; count() == 1 confirms the large bit is gone.
    CHECK(target.count() == 1);
}
