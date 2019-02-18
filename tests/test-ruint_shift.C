/* test_shift.cpp - Bits shift of RecInt test file

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
    ruint<STD_RECINT_SIZE> x, y, z;
    mpz_class size, gx, gy, gz, gcmp;
    USItype r;

    // Init. size = 2 ^ (2 ^ STD_RECINT_SIZE)
    mpz_ui_pow_ui(size.get_mpz_t(), 2, STD_RECINT_SIZE);
    mpz_ui_pow_ui(size.get_mpz_t(), 2, size.get_ui());
    RecInt::srand(limb(time(NULL)));

    // Loop
    for (UDItype l = 1; l < LOOPS; l++) {
        // Left shift
        rand(x);
        r = USItype(rand()) % (2 * NBBITS<STD_RECINT_SIZE>::value);
        ruint_to_mpz(gx, x);
        y = x << (UDItype)r; gy = gx << r; gy %= size;
        ruint_to_mpz(gcmp, y);
        if (gcmp != gy) return 1;

        x <<= (UDItype)r; gx <<= r; gx %= size;
        ruint_to_mpz(gcmp, x);
        if (gcmp != gx) return 1;

        // Right shift
        rand(x);
        r = USItype(rand()) % (2 * NBBITS<STD_RECINT_SIZE>::value);
        ruint_to_mpz(gx, x);
        y = x >> (UDItype)r; gy = gx >> r;
        ruint_to_mpz(gcmp, y);
        if (gcmp != gy) return 2;

        x >>= (UDItype)r; gx >>= r;
        ruint_to_mpz(gcmp, x);
        if (gcmp != gx) return 2;
    }

    return 0;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
