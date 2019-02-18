/* test_arith.cpp - Arithmetic operators of RecInt generic test file

   Return value.
   0    No error
   != 0 Bad result for an operation

   The following constants have to be defined.
   STD_RECINT_SIZE     size of recint (> 5)
   LOOPS           number of loops of randized tests
   RI_OP           RecInt arithmetic operation (RI_add_nc, RI_sub_nc, RI_mul)
   GMP_OP          GMP arithmetic operation (mpz_add, mpz_sub, mpz_mul)
   NOT_IN_PLACE    (optional) the in-place tests do not occur
   */

#include <cstddef> // required by gmp versions <= 5.1.3
#include <gmpxx.h>
#include <recint/recint.h>

#if not defined(LOOPS)
#define LOOPS 10000
#endif

using namespace RecInt;

int main(void)
{
    rmint<STD_RECINT_SIZE> x, y, z;
    mpz_class size, gx, gy, gz, gxy;
    USItype r;

    // Init. size = p
    RecInt::srand(limb(time(NULL)));
    ruint<STD_RECINT_SIZE> p;

    // Loop
    for (UDItype l = 1; l < LOOPS; l++) {
        // Init module
        do { rand(p); } while(p % 2 == 0);
        rmint<STD_RECINT_SIZE>::init_module(p);
        ruint_to_mpz(size, p);

        // Random sur rmint
        rand(x);
        rand(y);
        rand(z);

        // rmint_to_mpz x, y and the result to GMP
        rmint_to_mpz(gx, x);
        rmint_to_mpz(gy, y);
        rmint_to_mpz(gxy, z); // For addmul to be OK

        // Do RecInt operation
        RI_OP(z, x, y);
        rmint_to_mpz(gz, z);

        // Do GMP operation
        GMP_OP(gxy.get_mpz_t(), gx.get_mpz_t(), gy.get_mpz_t());
        gxy %= size; if (gxy < 0) gxy += size;

        // Compare both results
        if (gz != gxy) return 1;


        // Second test: with unsigned int
        r = USItype(rand());
        rmint_to_mpz(gx, x);

        RI_OP(z, x, r);
        rmint_to_mpz(gz, z);

        GMP_OPUI(gxy.get_mpz_t(), gx.get_mpz_t(), r);
        gxy %= size; if (gxy < 0) gxy += size;

        if (gz != gxy) return 2;


        // Third test: with repetition
        rmint_to_mpz(gx, x);
        RI_OP(x, x, x);

        rmint_to_mpz(gz, x);
        GMP_OP(gx.get_mpz_t(), gx.get_mpz_t(), gx.get_mpz_t());

        gx %= size; if (gx < 0) gx += size;
        if (gz != gx) return 3;


#if not defined(NOT_IN_PLACE)
        // Fourth test: in place
        rand(x);
        rand(y);
        rmint_to_mpz(gx, x);
        rmint_to_mpz(gy, y);

        RI_OP(x, y);
        rmint_to_mpz(gz, x);
        GMP_OP(gxy.get_mpz_t(), gx.get_mpz_t(), gy.get_mpz_t());

        gxy %= size; if (gxy < 0) gxy += size;
        if (gz != gxy) return 4;


        // Fifth test: in place with unsigned int
        r = USItype(rand());
        rmint_to_mpz(gx, x);

        RI_OP(x, r);
        rmint_to_mpz(gz, x);

        GMP_OPUI(gxy.get_mpz_t(), gx.get_mpz_t(), r);
        gxy %= size; if (gxy < 0) gxy += size;

        if (gz != gxy) return 5;


        // Sixth test: in place with repetition
        rmint_to_mpz(gx, x);
        RI_OP(x, x);

        rmint_to_mpz(gz, x);
        GMP_OP(gxy.get_mpz_t(), gx.get_mpz_t(), gx.get_mpz_t());

        gxy %= size; if (gxy < 0) gxy += size;
        if (gz != gxy) return 6;
#endif
    }

    return 0;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
