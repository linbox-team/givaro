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
#include <recint/ruint.h>

#if not defined(LOOPS)
#define LOOPS 10000
#endif

using namespace RecInt;

int main(void)
{
    ruint<STD_RECINT_SIZE> x, y, z;
    mpz_class size, gx, gy, gz, gxy;
    USItype r;

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

        // ruint_to_mpz x, y and the result to GMP
        ruint_to_mpz(gx, x);
        ruint_to_mpz(gy, y);
        ruint_to_mpz(gxy, z);

        // Do RecInt operation
        if (x > y) RI_OP(z, x, y);
        else RI_OP(z, y, x);
        ruint_to_mpz(gz, z);

        // Do GMP operation
        if (gx > gy) GMP_OP(gxy.get_mpz_t(), gx.get_mpz_t(), gy.get_mpz_t());
        else GMP_OP(gxy.get_mpz_t(), gy.get_mpz_t(), gx.get_mpz_t());
        gxy %= size;

        // Compare both results
        if (gz != gxy) return 1;


        // Second test: with unsigned int
        r = USItype(rand());
        while (x < r) { rand(x); r = USItype(rand()); }
        RI_OP(z, x, r);
        ruint_to_mpz(gx, x);
        ruint_to_mpz(gz, z);

        GMP_OPUI(gxy.get_mpz_t(), gx.get_mpz_t(), r);
        gxy %= size;
        if (gz != gxy) return 2;


        // Third test: with repetition
        while (x < 2) rand(x);

        ruint_to_mpz(gx, x);
        RI_OP(x, x, x);

        ruint_to_mpz(gz, x);
        GMP_OP(gx.get_mpz_t(), gx.get_mpz_t(), gx.get_mpz_t());
        gx %= size;
        if (gz != gx) return 3;

#if not defined(NOT_IN_PLACE)
        // Fourth test: in place
        rand(x);
        rand(y);
        ruint_to_mpz(gx, x);
        ruint_to_mpz(gy, y);

        if (x > y) {
            RI_OP(x, y);
            ruint_to_mpz(gz, x);
            GMP_OP(gxy.get_mpz_t(), gx.get_mpz_t(), gy.get_mpz_t());
        } else {
            RI_OP(y, x);
            ruint_to_mpz(gz, y);
            GMP_OP(gxy.get_mpz_t(), gy.get_mpz_t(), gx.get_mpz_t());
        }

        gxy %= size;
        if (gz != gxy) return 4;


        // Fifth test: in place with unsigned int
        r = USItype(rand());
        while (x < r) { rand(x); r = USItype(rand()); }
        ruint_to_mpz(gx, x);

        RI_OP(x, r);
        ruint_to_mpz(gz, x);

        GMP_OPUI(gxy.get_mpz_t(), gx.get_mpz_t(), r);
        gxy %= size;
        if (gz != gxy) return 5;


        // Sixth test: in place with repetition
        while (x < 2) rand(x);

        ruint_to_mpz(gx, x);
        RI_OP(x, x);

        ruint_to_mpz(gz, x);
        GMP_OP(gxy.get_mpz_t(), gx.get_mpz_t(), gx.get_mpz_t());
        gxy %= size;
        if (gz != gxy) return 6;
#endif
    }

    return 0;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
