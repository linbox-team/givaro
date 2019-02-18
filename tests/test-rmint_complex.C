/* recint_complex.cpp - Combined calculus with operators

   Return value.
   0    No error
   != 0 Bad result for an operation

   The following constants have to be defined.
   STD_RECINT_SIZE     gp of recint (> 5)
   LOOPS           number of loops of randomized tests
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
    ruint<STD_RECINT_SIZE> p;
    rmint<STD_RECINT_SIZE> a, b, c, d, z;
    mpz_class ga, gb, gc, gd, gz, gp, gcmp;

    // Init.
    RecInt::srand(limb(time(NULL)));
    do { rand(p); } while (p % 2 == 0);
    a.init_module(p);
    ruint_to_mpz(gp, p);

    // Loop
    for (UDItype l = 1; l < LOOPS; l++) {
        rand(a); rand(b); rand(c); rand(d);
        rmint_to_mpz(ga, a); rmint_to_mpz(gb, b);
        rmint_to_mpz(gc, c); rmint_to_mpz(gd, d);

        //--------- Tests rmint ----------

        // Test complex 1
        z = a * (b + c * d);
        gz = ga * (gb + gc * gd); gz %= gp;
        rmint_to_mpz(gcmp, z);
        if (gcmp != gz) return 1;

        // Test complex 2 (in place)
        a += a * (b + a);
        ga += ga * (gb + ga); ga %= gp;
        rmint_to_mpz(gcmp, a);
        if (gcmp != ga) return 2;

        // Test int * +
        z = 4 * a + b * (c + 1);
        gz = 4 * ga + gb * (gc + 1); gz %= gp;
        rmint_to_mpz(gcmp, z);
        if (gcmp != gz) return 3;

        // Test int -
        if (a > b) {
            z = c * (a - (b - 1));
            gz = gc * (ga - (gb - 1));
        } else {
            z = c * (b - (a - 1));
            gz = gc * (gb - (ga - 1));
        }
        gz %= gp; rmint_to_mpz(gcmp, z);
        if (gcmp != gz) return 4;

        // Tests in place et int
        a += 5; ga += 5; rmint_to_mpz(gcmp, a);
        if (a > 5 && gcmp != ga) return 5;

        if (b > 5) { b -= 5; gb -= 5; rmint_to_mpz(gcmp, b);
            if (gcmp != gb) return 6; }

        c *= 5; gc *= 5; gc %= gp; rmint_to_mpz(gcmp, c);
        if (gcmp != gc) return 7;
    }

    return 0;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
