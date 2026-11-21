#include <gss/innards/svo_bitset.hh>

#include <catch2/catch_test_macros.hpp>

#include <iterator>
#include <random>
#include <ranges>
#include <set>
#include <vector>

using namespace gss::innards;

namespace
{
    // SVOBitset stores up to 16 64-bit words inline ("short") and spills to the heap
    // ("long") beyond that. Every behavioural test is run at a small size (<= 1024
    // bits, inline) and a large size (> 1024 bits, heap) so both representations and
    // the boundary between them are exercised.
    constexpr unsigned small_bits = 128; // 2 words, inline
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

TEST_CASE("SVOBitset copy-assignment between heap bitsets of different sizes")
{
    // Assigning a larger heap bitset to a smaller one must reallocate to the larger size.
    // It used to reallocate to the old (smaller) size and then copy the larger one's words,
    // overflowing the buffer (caught by the sanitize build).
    SVOBitset larger(4 * large_bits, 0);
    larger.set(4 * large_bits - 1);
    SVOBitset smaller(large_bits, 0);
    smaller.set(3);

    smaller = larger;
    CHECK(smaller.test(4 * large_bits - 1));
    CHECK_FALSE(smaller.test(3));
    CHECK(smaller.count() == 1);

    SVOBitset back(large_bits, 0);
    back.set(7);
    smaller = back;
    CHECK(smaller.test(7));
    CHECK(smaller.count() == 1);
}

static_assert(std::forward_iterator<SVOBitset::Iterator>);
static_assert(std::ranges::forward_range<SVOBitset>);
static_assert(std::ranges::forward_range<const SVOBitset>);

namespace
{
    // Sizes either side of every boundary the representation has: within a word, at a word
    // edge, at the inline/heap switch (16 words = 1024 bits), and well onto the heap.
    const std::vector<unsigned> boundary_sizes = {1, 2, 63, 64, 65, 127, 128, 129, 1000, 1023, 1024, 1025, 1088, 2048, 3001};

    auto members(const SVOBitset & b) -> std::vector<unsigned>
    {
        std::vector<unsigned> result;
        for (auto v : b)
            result.push_back(v);
        return result;
    }

    auto random_set(unsigned size, double density, std::mt19937 & rand) -> std::pair<SVOBitset, std::set<unsigned>>
    {
        SVOBitset b(size, 0);
        std::set<unsigned> reference;
        std::bernoulli_distribution coin(density);
        for (unsigned v = 0; v < size; ++v)
            if (coin(rand)) {
                b.set(v);
                reference.insert(v);
            }
        return {b, reference};
    }
}

TEST_CASE("SVOBitset iteration of empty sets")
{
    SVOBitset default_constructed;
    CHECK(default_constructed.begin() == default_constructed.end());
    CHECK(default_constructed.find_last() == SVOBitset::npos);

    for (auto size : boundary_sizes) {
        SVOBitset b(size, 0);
        CHECK(b.begin() == b.end());
        CHECK(members(b).empty());
        CHECK(b.find_last() == SVOBitset::npos);
    }
}

TEST_CASE("SVOBitset iteration at word and storage boundaries")
{
    for (auto size : boundary_sizes) {
        SVOBitset b(size, 0);
        std::vector<unsigned> expected;
        for (unsigned v : {0u, 1u, 62u, 63u, 64u, 65u, 127u, 128u, 1023u, 1024u, size - 1})
            if (v < size && (expected.empty() || expected.back() < v)) {
                b.set(v);
                expected.push_back(v);
            }
        CHECK(members(b) == expected);
        CHECK(b.find_first() == expected.front());
        CHECK(b.find_last() == expected.back());
    }
}

TEST_CASE("SVOBitset iteration matches a reference set")
{
    std::mt19937 rand(2026);
    for (auto size : boundary_sizes)
        for (double density : {0.01, 0.1, 0.5, 0.9, 1.0}) {
            auto [b, reference] = random_set(size, density, rand);
            CHECK(members(b) == std::vector<unsigned>(reference.begin(), reference.end()));
            CHECK(std::vector<unsigned>(b.begin(), b.end()) == std::vector<unsigned>(reference.begin(), reference.end()));
            CHECK(std::ranges::distance(b) == std::ptrdiff_t(b.count()));
            CHECK(b.find_first() == (reference.empty() ? SVOBitset::npos : *reference.begin()));
            CHECK(b.find_last() == (reference.empty() ? SVOBitset::npos : *reference.rbegin()));
        }
}

TEST_CASE("SVOBitset iterators are multi-pass")
{
    SVOBitset b(large_bits, 0);
    b.set(3);
    b.set(700);
    b.set(1500);
    auto i = b.begin();
    auto j = i;
    CHECK(*i++ == 3);
    CHECK(*i == 700);
    CHECK(*j == 3); // copies advance independently
    CHECK(++j == i);
    CHECK(std::ranges::next(b.begin(), 3) == b.end());
}

TEST_CASE("SVOBitset reset_if removes what the predicate selects, in order")
{
    std::mt19937 rand(1854);
    for (auto size : boundary_sizes)
        for (double density : {0.05, 0.5, 1.0}) {
            auto [b, reference] = random_set(size, density, rand);

            std::vector<unsigned> visited;
            auto removed = b.reset_if([&](unsigned v) {
                visited.push_back(v);
                return v % 3 == 0 || v % 64 == 63;
            });

            // called once per member, in increasing order
            CHECK(visited == std::vector<unsigned>(reference.begin(), reference.end()));

            unsigned expected_removed = 0;
            std::vector<unsigned> expected_kept;
            for (auto v : reference)
                if (v % 3 == 0 || v % 64 == 63)
                    ++expected_removed;
                else
                    expected_kept.push_back(v);
            CHECK(removed == expected_removed);
            CHECK(members(b) == expected_kept);
            CHECK(b.count() == expected_kept.size());
        }
}

TEST_CASE("SVOBitset reset_below and reset_above")
{
    std::mt19937 rand(1066);
    for (auto size : boundary_sizes) {
        auto [full, reference] = random_set(size, 0.7, rand);
        for (unsigned a : {0u, 1u, 63u, 64u, 65u, size / 2, size - 1, size, size + 1, 5000u}) {
            SVOBitset below = full;
            below.reset_below(a);
            std::vector<unsigned> expected_below;
            for (auto v : reference)
                if (v >= a)
                    expected_below.push_back(v);
            CHECK(members(below) == expected_below);
            CHECK(below.count() == expected_below.size());

            SVOBitset above = full;
            above.reset_above(a);
            std::vector<unsigned> expected_above;
            for (auto v : reference)
                if (v <= a)
                    expected_above.push_back(v);
            CHECK(members(above) == expected_above);
            CHECK(above.count() == expected_above.size());
        }
    }
}
