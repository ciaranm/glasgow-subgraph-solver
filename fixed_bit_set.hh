/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_BIT_FIXED_BIT_SET_HH
#define GLASGOW_SUBGRAPH_SOLVER_BIT_FIXED_BIT_SET_HH 1

#include <array>
#include <vector>
#include <tuple>
#include <utility>
#include <algorithm>
#include <limits>

/// We'll use an array of unsigned long longs to represent our bits.
using BitWord = unsigned long long;

/// Number of bits per word.
static const constexpr int bits_per_word = sizeof(BitWord) * 8;

/**
 * A bitset with a fixed maximum size. This only provides the operations
 * we actually use in the bitset algorithms: it's more readable this way
 * than doing all the bit voodoo inline.
 *
 * Indices start at 0.
 */
template <unsigned words_>
class FixedBitSet
{
    private:
        using Bits = std::array<BitWord, words_>;

        Bits _bits = {{ }};

    public:
        static constexpr const unsigned npos = std::numeric_limits<unsigned>::max();

        FixedBitSet(unsigned, unsigned)
        {
        }

        FixedBitSet() = default;

        FixedBitSet(const FixedBitSet &) = default;

        FixedBitSet(FixedBitSet &&) = default;

        FixedBitSet & operator= (const FixedBitSet &) = default;

        FixedBitSet & operator= (FixedBitSet &&) = default;

        /**
         * Set a given bit 'on'.
         */
        auto set(int a) -> void
        {
            // The 1 does have to be of type BitWord. If we just specify a
            // literal, it ends up being an int, and it isn't converted
            // upwards until after the shift is done.
            _bits[a / bits_per_word] |= (BitWord{ 1 } << (a % bits_per_word));
        }

        /**
         * Set a given bit 'off'.
         */
        auto reset(int a) -> void
        {
            _bits[a / bits_per_word] &= ~(BitWord{ 1 } << (a % bits_per_word));
        }

        /**
         * Set all bits off.
         */
        auto reset() -> void
        {
            for (unsigned i = 0 ; i < words_ ; ++i)
                _bits[i] = 0;
        }

        /**
         * Is a given bit on?
         */
        auto test(int a) const -> bool
        {
            return _bits[a / bits_per_word] & (BitWord{ 1 } << (a % bits_per_word));
        }

        /**
         * How many bits are on?
         */
        auto count() const -> unsigned
        {
            unsigned result = 0;
            for (auto & p : _bits)
                result += __builtin_popcountll(p);
            return result;
        }

        /**
         * Intersect (bitwise-and) with another set.
         */
        auto operator &= (const FixedBitSet<words_> & other) -> void
        {
            for (typename Bits::size_type i = 0 ; i < words_ ; ++i)
                _bits[i] = _bits[i] & other._bits[i];
        }

        /**
         * Union (bitwise-or) with another set.
         */
        auto operator|= (const FixedBitSet<words_> & other) -> void
        {
            for (typename Bits::size_type i = 0 ; i < words_ ; ++i)
                _bits[i] = _bits[i] | other._bits[i];
        }

        auto operator~ () -> FixedBitSet<words_> const
        {
            FixedBitSet<words_> result;
            for (typename Bits::size_type i = 0 ; i < words_ ; ++i)
                result._bits[i] = ~_bits[i];
            return result;
        }

        /**
         * Return the index of the first set ('on') bit, or npos if we are
         * empty.
         */
        auto find_first() const -> unsigned
        {
            for (typename Bits::size_type i = 0 ; i < _bits.size() ; ++i) {
                int b = __builtin_ffsll(_bits[i]);
                if (0 != b)
                    return i * bits_per_word + b - 1;
            }
            return npos;
        }
};

#endif
