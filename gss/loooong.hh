#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_LOOOONG_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_LOOOONG_HH 1

#include <cstring>
#include <gmp.h>
#include <iostream>

using std::strlen;

namespace gss
{
    class loooong
    {
    private:
        mpz_t value;

    public:
        loooong()
        {
            mpz_init(value);
        }
        loooong(long v)
        {
            mpz_init_set_si(value, v);
        }
        loooong(unsigned long v)
        {
            mpz_init_set_ui(value, v);
        }
        loooong(int v)
        {
            mpz_init_set_si(value, v);
        }
        loooong(const loooong & other)
        {
            mpz_init_set(value, other.value);
        }
        loooong(loooong && other) noexcept
        {
            mpz_init(value);
            mpz_set(value, other.value);
        }

        ~loooong()
        {
            mpz_clear(value);
        }

        mpz_t & raw()
        {
            return value;
        }
        const mpz_t & raw() const
        {
            return value;
        }

        loooong & operator=(const loooong & other)
        {
            if (this != &other) mpz_set(value, other.value);
            return *this;
        }

        loooong & operator=(long v)
        {
            mpz_set_si(value, v);
            return *this;
        }

        loooong & operator=(unsigned long v)
        {
            mpz_set_ui(value, v);
            return *this;
        }

        loooong operator+(const loooong & rhs) const
        {
            loooong r;
            mpz_add(r.value, value, rhs.value);
            return r;
        }

        loooong operator-(const loooong & rhs) const
        {
            loooong r;
            mpz_sub(r.value, value, rhs.value);
            return r;
        }

        loooong operator*(const loooong & rhs) const
        {
            loooong r;
            mpz_mul(r.value, value, rhs.value);
            return r;
        }

        loooong & operator+=(const loooong & rhs)
        {
            mpz_add(value, value, rhs.value);
            return *this;
        }

        loooong & operator-=(const loooong & rhs)
        {
            mpz_sub(value, value, rhs.value);
            return *this;
        }

        loooong & operator*=(const loooong & rhs)
        {
            mpz_mul(value, value, rhs.value);
            return *this;
        }

        loooong & operator++()
        {
            mpz_add_ui(value, value, 1);
            return *this;
        }

        loooong & operator/=(const loooong & rhs)
        {
            mpz_tdiv_q(value, value, rhs.value);
            return *this;
        }

        loooong operator/(const loooong & rhs) const
        {
            loooong r;
            mpz_tdiv_q(r.raw(), value, rhs.value);
            return r;
        }

        std::strong_ordering operator<=>(const loooong& rhs) const noexcept {
            int cmp = mpz_cmp(value, rhs.value);
            if (cmp < 0) return std::strong_ordering::less;
            if (cmp > 0) return std::strong_ordering::greater;
            return std::strong_ordering::equal;
        }

        friend std::ostream & operator<<(std::ostream & os, const loooong & x)
        {
            char * str = mpz_get_str(nullptr, 10, x.value);
            os << str;
            void (*freefunc)(void *, size_t);
            mp_get_memory_functions(nullptr, nullptr, &freefunc);
            freefunc(str, std::strlen(str) + 1);
            return os;
        }
    };

    inline loooong gcd(const loooong & a, const loooong & b)
    {
        loooong r;
        mpz_gcd(r.raw(), a.raw(), b.raw());
        return r;
    }

}

#endif
