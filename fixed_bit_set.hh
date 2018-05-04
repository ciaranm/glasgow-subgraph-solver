/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GUARD_BIT_GRAPH_HH
#define GUARD_BIT_GRAPH_HH 1

#include <array>
#include <vector>
#include <tuple>
#include <utility>
#include <algorithm>

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
        auto unset(int a) -> void
        {
            _bits[a / bits_per_word] &= ~(BitWord{ 1 } << (a % bits_per_word));
        }

        /**
         * Set all bits off.
         */
        auto unset_all() -> void
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
        auto popcount() const -> unsigned
        {
            unsigned result = 0;
            for (auto & p : _bits)
                result += __builtin_popcountll(p);
            return result;
        }

        /**
         * Intersect (bitwise-and) with another set.
         */
        auto intersect_with(const FixedBitSet<words_> & other) -> void
        {
            for (typename Bits::size_type i = 0 ; i < words_ ; ++i)
                _bits[i] = _bits[i] & other._bits[i];
        }

        /**
         * Union (bitwise-or) with another set.
         */
        auto union_with(const FixedBitSet<words_> & other) -> void
        {
            for (typename Bits::size_type i = 0 ; i < words_ ; ++i)
                _bits[i] = _bits[i] | other._bits[i];
        }

        /**
         * Intersect with the complement of another set.
         */
        auto intersect_with_complement(const FixedBitSet<words_> & other) -> void
        {
            for (typename Bits::size_type i = 0 ; i < words_ ; ++i)
                _bits[i] = _bits[i] & ~other._bits[i];
        }

        /**
         * Return the index of the first set ('on') bit, or -1 if we are
         * empty.
         */
        auto first_set_bit() const -> int
        {
            for (typename Bits::size_type i = 0 ; i < _bits.size() ; ++i) {
                int b = __builtin_ffsll(_bits[i]);
                if (0 != b)
                    return i * bits_per_word + b - 1;
            }
            return -1;
        }

        /**
         * Return the index of the first set ('on') bit, or -1 if we are
         * empty. Start at word offset, which is updated.
         */
        auto first_set_bit_from(unsigned & offset) const -> int
        {
            for ( ; offset < _bits.size() ; ++offset) {
                int b = __builtin_ffsll(_bits[offset]);
                if (0 != b)
                    return offset * bits_per_word + b - 1;
            }
            return -1;
        }
};

/**
 * We have to decide at compile time what the largest graph we'll support
 * is.
 */
constexpr auto max_graph_words __attribute__((unused)) = 1024;

/**
 * Thrown if we exceed max_graph_words.
 */
class GraphTooBig :
    public std::exception
{
    public:
        GraphTooBig() throw ();
};

#endif
