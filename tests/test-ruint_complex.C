/* recint_complex.cpp - Combined calculus with operators

   Return value.
   0    No error
   != 0 Bad result for an operation

   The following constants have to be defined.
   STD_RECINT_SIZE     size of recint (> 5)
   LOOPS           number of loops of randomized tests
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
    ruint<STD_RECINT_SIZE+1> zl;
    ruint<STD_RECINT_SIZE> a, b, c, d, z;
    mpz_class size, ga, gb, gc, gd, gz, gcmp;

    // Init. size = 2 ^ (2 ^ STD_RECINT_SIZE)
    mpz_ui_pow_ui(size.get_mpz_t(), 2, STD_RECINT_SIZE);
    mpz_ui_pow_ui(size.get_mpz_t(), 2, size.get_ui());
    RecInt::srand(limb(time(NULL)));

    a = 0;
    b = -a;

    if (a != b)
        return -1;

    // Loop
    for (UDItype l = 1; l < LOOPS; l++) {
        rand(a); rand(b); rand(c); rand(d);
        ruint_to_mpz(ga, a); ruint_to_mpz(gb, b);
        ruint_to_mpz(gc, c); ruint_to_mpz(gd, d);

        //--------- Tests ruint ----------

        // Test complex 1
        z = a * (b + c * d);
        gz = ga * (gb + gc * gd); gz %= size;
        ruint_to_mpz(gcmp, z);
        if (gcmp != gz) return 1;

        // Test complex 2 (in place)
        a += a * (b + a);
        ga += ga * (gb + ga); ga %= size;
        ruint_to_mpz(gcmp, a);
        if (gcmp != ga) return 2;

        // Test complex 3 (shift)
        z = (a * z) << (2 * STD_RECINT_SIZE);
        gz = (ga * gz) << (2 * STD_RECINT_SIZE); gz %= size;
        ruint_to_mpz(gcmp, z);
        if (gcmp != gz) return 3;

        // Test int * +
        z = 4 * a + b * (c + 1);
        gz = 4 * ga + gb * (gc + 1); gz %= size;
        ruint_to_mpz(gcmp, z);
        if (gcmp != gz) return 4;

        // Test int -
        if (a > b) {
            z = c * (a - (b - 1));
            gz = gc * (ga - (gb - 1));
        } else {
            z = c * (b - (a - 1));
            gz = gc * (gb - (ga - 1));
        }
        gz %= size; ruint_to_mpz(gcmp, z);
        if (gcmp != gz) return 5;

        // Test int / %
        z = b + (a / 2) % 3;
        gz = gb + (ga / 2) % 3;
        ruint_to_mpz(gcmp, z);
        if (gcmp != gz) return 6;

        // Tests in place et int
        a += 5; ga += 5; ruint_to_mpz(gcmp, a);
        if (a > 5 && gcmp != ga) return 7;

        if (b > 5) { b -= 5; gb -= 5; ruint_to_mpz(gcmp, b);
            if (gcmp != gb) return 8; }

        c *= 5; gc *= 5; gc %= size; ruint_to_mpz(gcmp, c);
        if (gcmp != gc) return 9;

        d /= 5; gd /= 5; ruint_to_mpz(gcmp, d);
        if (gcmp != gd) return 10;

        a %= 5; ga %= 5; ruint_to_mpz(gcmp, a);
        if (gcmp != ga) return 11;
    }

    return 0;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
