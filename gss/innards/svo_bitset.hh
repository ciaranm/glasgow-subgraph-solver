#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_SVO_BITSET_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_SVO_BITSET_HH 1

#include <algorithm>
#include <array>
#include <cstring>
#include <limits>

#ifndef USE_PORTABLE_SNIPPETS_BUILTIN
#include <bit>
#endif

namespace gss::innards
{
#ifdef USE_PORTABLE_SNIPPETS_BUILTIN
#include <portable-snippets/builtin/builtin.h>
    int popcount(unsigned long long x)
    {
        return psnip_builtin_popcountll(x);
    }
    int countr_zero(unsigned long long x)
    {
        return psnip_builtin_ctzll(x);
    }
#else
    using std::countr_zero;
    using std::popcount;
#endif

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

    public:
        static constexpr const unsigned npos = std::numeric_limits<unsigned>::max();

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

        auto find_first() const -> unsigned
        {
            const BitWord * b = (_is_long() ? _data.long_data : _data.short_data);
            for (unsigned i = 0; i < n_words; ++i) {
                if (0 != b[i]) // test the input to countr_zero because on some systems, the output is undefined at b[i]=0
                    return i * bits_per_word + countr_zero(b[i]);
            }
            return npos;
        }

        auto reset(int a) -> void
        {
            BitWord * b = (_is_long() ? _data.long_data : _data.short_data);
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
            BitWord * b = (_is_long() ? _data.long_data : _data.short_data);
            b[a / bits_per_word] |= (BitWord{1} << (a % bits_per_word));
        }

        auto test(int a) const -> bool
        {
            const BitWord * b = (_is_long() ? _data.long_data : _data.short_data);
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
            const BitWord * b = (_is_long() ? _data.long_data : _data.short_data);
            for (unsigned i = 0, i_end = n_words; i < i_end; ++i)
                result += popcount(b[i]);

            return result;
        }
    };
}

#endif
