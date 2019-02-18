/* recint_extra.cpp - Extra arithmetics function of RecInt library test file

   Return value.
   0    No error
   != 0 Bad result for an operation

   The following constants have to be defined.
   STD_RECINT_SIZE     size of recint (> 5)
   LOOPS           number of loops of randomized tests
   */

#include <cstddef> // required by gmp versions <= 5.1.3
#include <gmpxx.h>
#include <recint/recint.h>

#if not defined(LOOPS)
#define LOOPS 1000
#endif

using namespace RecInt;

int main(void)
{
    ruint<STD_RECINT_SIZE> a, b, z, p;
    rmint<STD_RECINT_SIZE> ma, mb, mz;
    mpz_class ga, gb, gz, gcmp;

    // Init.
    RecInt::srand(limb(time(NULL)));

    // Loop
    for (UDItype l = 1; l < LOOPS; l++) {
        rand(a); rand(b);
        ruint_to_mpz(ga, a); ruint_to_mpz(gb, b);

        // GCD in ruint
        gcd(z, a, b);
        mpz_gcd(gz.get_mpz_t(), ga.get_mpz_t(), gb.get_mpz_t());
        ruint_to_mpz(gcmp, z);
        if (gcmp != gz) return 1;

        // rmint Modular square root (if quadratic residue)
        do { rand(p); } while(p % 2 == 0);
        ma.init_module(p);
        rand(ma);
        square_root(mz, ma);
        if (mz == 0 && is_quadratic_residue(ma)) return 2;
        else if (mz != 0 && !is_quadratic_residue(ma)) return 3;
        else if (mz != 0 && (mz * mz) != ma) return 4;
    }

    return 0;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
