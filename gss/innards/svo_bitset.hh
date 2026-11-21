#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_SVO_BITSET_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_SVO_BITSET_HH 1

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstring>
#include <iterator>
#include <limits>

#ifndef USE_PORTABLE_SNIPPETS_BUILTIN
#include <bit>
#endif

namespace gss::innards
{
#ifdef USE_PORTABLE_SNIPPETS_BUILTIN
#include <portable-snippets/builtin/builtin.h>
    inline int popcount(unsigned long long x)
    {
        return psnip_builtin_popcountll(x);
    }
    inline int countr_zero(unsigned long long x)
    {
        return psnip_builtin_ctzll(x);
    }
    inline int countl_zero(unsigned long long x)
    {
        return psnip_builtin_clzll(x);
    }
#else
    using std::countl_zero;
    using std::countr_zero;
    using std::popcount;
#endif

    // A fixed-size set of small unsigned integers, stored as a bitset with inline storage
    // for up to svo_size words.
    //
    // Iterating a set visits its members in increasing order:
    //
    //     for (auto v : set) ...
    //
    // The iterator is a std::forward_iterator, so the standard algorithms and range
    // constructors work too (vector<int> members(set.begin(), set.end())). Do not modify a
    // set while iterating over it; to remove members selected by a predicate, use reset_if,
    // which visits members in the same order.
    class SVOBitset
    {
    private:
        using BitWord = unsigned long long;
        static const constexpr int bits_per_word = sizeof(BitWord) * 8;
        static const constexpr int svo_size = 16;

        union {
            BitWord short_data[svo_size];
            BitWord * long_data;
        } _data;

        unsigned n_words;

        constexpr auto _is_long() const -> bool
        {
            return n_words > svo_size;
        }

        auto _words() -> BitWord *
        {
            return _is_long() ? _data.long_data : _data.short_data;
        }

        auto _words() const -> const BitWord *
        {
            return _is_long() ? _data.long_data : _data.short_data;
        }

    public:
        static constexpr const unsigned npos = std::numeric_limits<unsigned>::max();

        class Iterator
        {
        private:
            const BitWord * _words = nullptr;
            unsigned _n_words = 0;
            unsigned _word_index = 0;
            BitWord _remaining = 0; // the members of _words[_word_index] not yet visited

            // move to the next non-empty word, or to the end
            auto _advance() -> void
            {
                while (++_word_index < _n_words)
                    if (0 != (_remaining = _words[_word_index]))
                        return;
                _remaining = 0;
            }

            friend class SVOBitset;

            Iterator(const BitWord * words, unsigned n_words, unsigned word_index) :
                _words(words),
                _n_words(n_words),
                _word_index(word_index)
            {
                if (_word_index < _n_words && 0 == (_remaining = _words[_word_index]))
                    _advance();
            }

        public:
            using iterator_concept = std::forward_iterator_tag;
            using iterator_category = std::input_iterator_tag; // reference is not a real reference
            using value_type = unsigned;
            using difference_type = std::ptrdiff_t;
            using reference = unsigned;
            using pointer = void;

            Iterator() = default;

            auto operator*() const -> unsigned
            {
                return _word_index * bits_per_word + countr_zero(_remaining);
            }

            auto operator++() -> Iterator &
            {
                _remaining &= _remaining - 1;
                if (0 == _remaining)
                    _advance();
                return *this;
            }

            auto operator++(int) -> Iterator
            {
                auto result = *this;
                ++*this;
                return result;
            }

            // Only iterators over the same set are comparable, so the position is enough.
            auto operator==(const Iterator & other) const -> bool
            {
                return _word_index == other._word_index && _remaining == other._remaining;
            }
        };

        using iterator = Iterator;
        using const_iterator = Iterator;

        SVOBitset()
        {
            n_words = 0;
        }

        SVOBitset(unsigned size, unsigned bits);

        SVOBitset(const SVOBitset & other)
        {
            if (other._is_long()) {
                n_words = other.n_words;
                _data.long_data = new BitWord[n_words];
                std::copy(other._data.long_data, other._data.long_data + other.n_words, _data.long_data);
            }
            else {
                n_words = other.n_words;
                std::copy(&other._data.short_data[0], &other._data.short_data[svo_size], &_data.short_data[0]);
            }
        }

        ~SVOBitset()
        {
            if (_is_long())
                delete[] _data.long_data;
        }

        auto operator=(const SVOBitset & other) -> SVOBitset &
        {
            if (&other == this)
                return *this;

            if (other._is_long()) {
                if (! _is_long()) {
                    n_words = other.n_words;
                    _data.long_data = new BitWord[n_words];
                }
                else if (n_words != other.n_words) {
                    delete[] _data.long_data;
                    n_words = other.n_words;
                    _data.long_data = new BitWord[n_words];
                }

                std::copy(other._data.long_data, other._data.long_data + other.n_words, _data.long_data);
            }
            else {
                if (_is_long())
                    delete[] _data.long_data;
                n_words = other.n_words;
                std::copy(&other._data.short_data[0], &other._data.short_data[svo_size], &_data.short_data[0]);
            }

            return *this;
        }

        auto any() const -> bool
        {
            if (! _is_long()) {
                for (unsigned i = 0; i < n_words; ++i)
                    if (0 != _data.short_data[i])
                        return true;

                return false;
            }
            else {
                for (unsigned i = 0, i_end = n_words; i < i_end; ++i)
                    if (0 != _data.long_data[i])
                        return true;

                return false;
            }
        }

        auto begin() const -> Iterator
        {
            return Iterator{_words(), n_words, 0};
        }

        auto end() const -> Iterator
        {
            return Iterator{_words(), n_words, n_words};
        }

        // The smallest member, or npos if the set is empty.
        auto find_first() const -> unsigned
        {
            const BitWord * b = _words();
            for (unsigned i = 0; i < n_words; ++i) {
                if (0 != b[i]) // test the input to countr_zero because on some systems, the output is undefined at b[i]=0
                    return i * bits_per_word + countr_zero(b[i]);
            }
            return npos;
        }

        // The largest member, or npos if the set is empty.
        auto find_last() const -> unsigned
        {
            const BitWord * b = _words();
            for (unsigned i = n_words; i > 0; --i) {
                if (0 != b[i - 1])
                    return i * bits_per_word - 1 - countl_zero(b[i - 1]);
            }
            return npos;
        }

        // Remove every member v for which pred(v) is true, and return how many were removed.
        // pred is called exactly once for each member, in increasing order, so it may have
        // side effects (such as logging a proof step for each removal). It must not modify
        // this set.
        template <typename Pred_>
        auto reset_if(Pred_ && pred) -> unsigned
        {
            BitWord * b = _words();
            unsigned removed = 0;
            for (unsigned i = 0; i < n_words; ++i) {
                BitWord to_remove = 0;
                for (BitWord w = b[i]; 0 != w; w &= w - 1) {
                    auto bit = countr_zero(w);
                    if (pred(i * bits_per_word + bit))
                        to_remove |= BitWord{1} << bit;
                }
                b[i] &= ~to_remove;
                removed += popcount(to_remove);
            }
            return removed;
        }

        // Remove every member less than a.
        auto reset_below(unsigned a) -> void
        {
            BitWord * b = _words();
            unsigned full_words = std::min(a / bits_per_word, n_words);
            std::fill(b, b + full_words, 0);
            if (full_words < n_words && 0 != a % bits_per_word)
                b[full_words] &= ~BitWord{0} << (a % bits_per_word);
        }

        // Remove every member greater than a.
        auto reset_above(unsigned a) -> void
        {
            if (a >= n_words * bits_per_word)
                return;
            BitWord * b = _words();
            unsigned word = a / bits_per_word, bit = a % bits_per_word;
            if (bit + 1 < bits_per_word)
                b[word] &= ~(~BitWord{0} << (bit + 1));
            std::fill(b + word + 1, b + n_words, 0);
        }

        auto reset(int a) -> void
        {
            BitWord * b = _words();
            b[a / bits_per_word] &= ~(BitWord{1} << (a % bits_per_word));
        }

        auto reset() -> void
        {
            if (! _is_long())
                std::fill(&_data.short_data[0], &_data.short_data[svo_size], 0);
            else
                std::fill(_data.long_data, _data.long_data + n_words, 0);
        }

        auto set(int a) -> void
        {
            BitWord * b = _words();
            b[a / bits_per_word] |= (BitWord{1} << (a % bits_per_word));
        }

        auto test(int a) const -> bool
        {
            const BitWord * b = _words();
            return b[a / bits_per_word] & (BitWord{1} << (a % bits_per_word));
        }

        auto operator&=(const SVOBitset & other) -> SVOBitset &
        {
            if (! _is_long()) {
                for (unsigned i = 0; i < svo_size; ++i)
                    _data.short_data[i] &= other._data.short_data[i];
            }
            else {
                for (unsigned i = 0; i < n_words; ++i)
                    _data.long_data[i] &= other._data.long_data[i];
            }

            return *this;
        }

        auto operator|=(const SVOBitset & other) -> SVOBitset &
        {
            if (! _is_long()) {
                for (unsigned i = 0; i < svo_size; ++i)
                    _data.short_data[i] |= other._data.short_data[i];
            }
            else {
                for (unsigned i = 0; i < n_words; ++i)
                    _data.long_data[i] |= other._data.long_data[i];
            }

            return *this;
        }

        auto intersect_with_complement(const SVOBitset & other) -> void
        {
            if (! _is_long()) {
                for (unsigned i = 0; i < svo_size; ++i)
                    _data.short_data[i] &= ~other._data.short_data[i];
            }
            else {
                for (unsigned i = 0; i < n_words; ++i)
                    _data.long_data[i] &= ~other._data.long_data[i];
            }
        }

        auto count() const -> unsigned
        {
            unsigned result = 0;
            const BitWord * b = _words();
            for (unsigned i = 0, i_end = n_words; i < i_end; ++i)
                result += popcount(b[i]);

            return result;
        }
    };
}

#endif
