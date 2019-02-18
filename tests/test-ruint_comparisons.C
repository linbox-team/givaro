/* ruint_comparisons.cpp - Comparison functions for ruint of RecInt library test file

   Return value.
   0    No error
   != 0 Bad result for an operation

   The following constants have to be defined.
   STD_RECINT_SIZE     size of recint (> 5)
   LOOPS           number of loops of randized tests
   */

#include <cstddef> // required by gmp versions <= 5.1.3
#include <gmpxx.h>
#include <recint/ruint.h>

#if not defined(LOOPS)
#define LOOPS 10000
#endif

using namespace RecInt;

int main(void)
{
    ruint<STD_RECINT_SIZE> x, y, zero(0), one(1);
    mpz_class gx, gy, gcmp;
    USItype r;

    // Init.
    RecInt::srand(limb(time(NULL)));

    // Loop
    for (UDItype l = 1; l < LOOPS; l++) {
        // RI rand
        rand(x);
        rand(y);
        ruint_to_mpz(gx, x);
        ruint_to_mpz(gy, y);

        // Comp function
        if (cmp(x, x) != 0) return 1;
        if (cmp(x, y) > 0 && gx <= gy) return 1;
        if (cmp(x, y) < 0 && gx >= gy) return 1;

        r = USItype(rand());
        if (cmp(x, (UDItype)r) == 0 && gx != r) return 1;
        if (cmp(x, (UDItype)r) > 0 && gx <= r) return 1;
        if (cmp(x, (UDItype)r) < 0 && gx >= r) return 1;

        r = USItype(rand());
        if (cmp(x, (USItype)r) == 0 && gx != r) return 1;
        if (cmp(x, (USItype)r) > 0 && gx <= r) return 1;
        if (cmp(x, (USItype)r) < 0 && gx >= r) return 1;

        r = USItype(rand());
        if (cmp(x, (DItype)r) == 0 && gx != r) return 1;
        if (cmp(x, (DItype)r) > 0 && gx <= r) return 1;
        if (cmp(x, (DItype)r) < 0 && gx >= r) return 1;

        r = USItype(rand());
        if (cmp(x, (SItype)r) == 0 && gx != r) return 1;
        if (cmp(x, (SItype)r) > 0 && gx <= r) return 1;
        if (cmp(x, (SItype)r) < 0 && gx >= r) return 1;

        // Misc functions
        if (x != x) return 2;
        if (x != y && gx == gy) return 2;

        if (!(x == x)) return 3;
        if (x == y && gx != gy) return 3;

        if (x != 0 && gx == 0) return 4;
        if (zero != 0) return 4;

        if (x == 0 && gx != 0) return 5;
        if (!(zero == 0)) return 5;

        if (x != 1 && gx == 1) return 6;
        if (one != 1) return 6;

        if (x == 1 && gx != 1) return 7;
        if (!(one == 1)) return 7;

        // Comp functions
        if (x > x) return 8;
        if (x > y && gx <= gy) return 8;

        if (!(x >= x)) return 9;
        if (x >= y && gx < gy) return 9;

        if (x < x) return 10;
        if (x < y && gx >= gy) return 10;

        if (!(x <= x)) return 11;
        if (x <= y && gx > gy) return 11;
    }

    return 0;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
