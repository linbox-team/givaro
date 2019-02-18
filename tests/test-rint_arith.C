/* test_arith.cpp - Arithmetic operators of RecInt generic test file

   Return value.
   0    No error
   != 0 Bad result for an operation

   The following constants have to be defined.
   STD_RECINT_SIZE     size of recint (> 5)
   LOOPS           number of loops of randized tests
   RI_OP           RecInt arithmetic operation (RI_add_nc, RI_sub_nc, RI_mul)
   GMP_OP          GMP arithmetic operation (mpz_add, mpz_sub, mpz_mul)
   NOT_IN_PLACE    (optional) the in-place tests do not occur [for addmul and div]
   */

#include <cstddef> // required by gmp versions <= 5.1.3
#include <gmpxx.h>
#include <recint/rint.h>

#include <recint/ruint.h>
#include <iostream>

#if not defined(LOOPS)
#define LOOPS 10000
#endif

using namespace RecInt;

int main(void)
{
    RecInt::rint<STD_RECINT_SIZE> x, y, z;
    mpz_class size, gx, gy, gz, gxy;

    // Init. size = 2 ^ (2 ^ STD_RECINT_SIZE)
    mpz_ui_pow_ui(size.get_mpz_t(), 2, STD_RECINT_SIZE);
    mpz_ui_pow_ui(size.get_mpz_t(), 2, size.get_ui());
    RecInt::srand(limb(time(NULL)));

    // Loop
    for (UDItype l = 1; l < LOOPS; l++) {
        // RecInt rand
        rand(x);
        rand(y);
        rand(z);

#if defined(ARITH_MOD)
        // Used by div and mod to be sure that the divisor is positive
        if (y.isNegative()) y = -y;
#endif

        // rint_to_mpz x, y and the result to GMP
        rint_to_mpz(gx, x);
        rint_to_mpz(gy, y);
        rint_to_mpz(gxy, z);

        RI_OP(z, x, y);
        GMP_OP(gxy.get_mpz_t(), gx.get_mpz_t(), gy.get_mpz_t());
        rint_to_mpz(gz, z);

        // Overflow correction
#if not defined(ARITH_MOD)
        gxy %= size;
        if (gxy < - size / 2u) gxy += size;
        if (gxy >= size / 2u)  gxy -= size;
#endif

        // Compare both results
        if (gz != gxy) return 1;

#if defined(GMP_OPUI)
        // Second test: with unsigned int
        USItype r = USItype(rand());
        while (x < r) { rand(x); r = USItype(rand()); }
        RI_OP(z, x, r);
        rint_to_mpz(gx, x);
        rint_to_mpz(gz, z);

        GMP_OPUI(gxy.get_mpz_t(), gx.get_mpz_t(), r);

        // Overflow correction
        gxy %= size;
        if (gxy < - size / 2u) gxy += size;
        if (gxy >= size / 2u)  gxy -= size;

        if (gz != gxy) return 2;
#endif


        // Third test: with repetition
#if defined(ARITH_MOD)
        while (x < 0) rand(x);
#else
        while (x == 0) rand(x);
#endif
        rint_to_mpz(gx, x);

        RI_OP(x, x, x);
        rint_to_mpz(gz, x);

        GMP_OP(gx.get_mpz_t(), gx.get_mpz_t(), gx.get_mpz_t());

        // Overflow correction
#if not defined(ARITH_MOD)
        gx %= size;
        if (gx < - size / 2u) gx += size;
        if (gx >= size / 2u)  gx -= size;
#endif

        if (gz != gx) return 3;

#if not defined(NOT_IN_PLACE)
        // Fourth test: in place
        rand(x);
        rand(y);
        rint_to_mpz(gx, x);
        rint_to_mpz(gy, y);

        RI_OP(x, y);
        rint_to_mpz(gz, x);
        GMP_OP(gxy.get_mpz_t(), gx.get_mpz_t(), gy.get_mpz_t());

        // Overflow correction
        gxy %= size;
        if (gxy < - size / 2u) gxy += size;
        if (gxy >= size / 2u)  gxy -= size;

        if (gz != gxy) return 4;


        // Fifth test: in place with unsigned int
        r = USItype(rand());
        while (x < r) { rand(x); r = USItype(rand()); }
        rint_to_mpz(gx, x);

        RI_OP(x, r);
        rint_to_mpz(gz, x);

        GMP_OPUI(gxy.get_mpz_t(), gx.get_mpz_t(), r);

        // Overflow correction
        gxy %= size;
        if (gxy < - size / 2u) gxy += size;
        if (gxy >= size / 2u)  gxy -= size;

        if (gz != gxy) return 5;


        // Sixth test: in place with repetition
        while (x < 2) rand(x);

        rint_to_mpz(gx, x);
        RI_OP(x, x);

        rint_to_mpz(gz, x);
        GMP_OP(gxy.get_mpz_t(), gx.get_mpz_t(), gx.get_mpz_t());

        // Overflow correction
        gxy %= size;
        if (gxy < - size / 2u) gxy += size;
        if (gxy >= size / 2u)  gxy -= size;

        if (gz != gxy) return 6;
#endif
    }

    return 0;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
