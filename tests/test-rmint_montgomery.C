/* rmint_montgomery.cpp - Montgomery algorithms of RecInt library test file

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
#define LOOPS 1000
#endif

// This file sets MG_DEFAULT to MG_ACTIVE
// but this is overwritten by the Makefile... so:
#undef MG_DEFAULT
#include <recint/rmintmg.h>

using namespace RecInt;

int main(void)
{
    rmint<STD_RECINT_SIZE> a, b, c, d;
    ruint<STD_RECINT_SIZE + 1> cl;
    ruint<STD_RECINT_SIZE> power;
    mpz_class ga, gb, gc, gp, gcmp;

    // Init. size = p
    RecInt::srand(limb(time(NULL)));
    ruint<STD_RECINT_SIZE> p;

    // Loop
    for (UDItype l = 1; l < LOOPS; l++) {
        do { rand(p); } while ((p % 2) == 0);
        a.init_module(p);
        ruint_to_mpz(gp, p);

        //------- Using direct functions --------//

        rand(b); // b must be in MG space
        reduction(a, b); // a should be such that 0 <= a < p
        if (a.Value >= p) return 1;

        to_mg(c, a);
        if (c != b) return 2;

        lmul(cl, a.Value, a.Value); // In classic space
        mod_n(c.Value, cl, c.p);
        mul(d, b, b); // In Montgomery space
        reduction(d);
        if (c != d) return 3;

        rand(power);
        exp_mod(c.Value, a.Value, power, c.p); // In classic space
        exp(d, b, power); // In Montgomery space
        reduction(d, d);
        if (c != d) return 4;

        //------- Using indirect functions --------//

        rand(a); rand(b); // a and b must be in MG space
        rmint_to_mpz(ga, a); rmint_to_mpz(gb, b); // ga and gb should have been reduced
        c = a * b; // c must be in MG space
        gc = ga * gb; gc %= gp; // Classic multiplication
        rmint_to_mpz(gcmp, c); // c should be the correct answer
        if (gcmp != gc) return 5;

        c = a + b; // c must be in MG space
        gc = ga + gb; gc %= gp; // Classic operation
        rmint_to_mpz(gcmp, c); // c should be the correct answer
        if (gcmp != gc) return 6;

        while (gcd(get_ruint(b), p) != 1) rand(b); // For b to be invertible
        c = a / b; // c must be in MG space
        c *= b;
        if (c != a) return 7;
    }

    return 0;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
