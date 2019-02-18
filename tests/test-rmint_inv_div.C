/* rmint_inv.cpp - Inverse in modular arithmetic of RecInt library test file

   Return value.
   0    No error
   != 0 Bad result for an operation

   The following constants have to be defined.
   STD_RECINT_SIZE     size of recint (> 5)
   LOOPS           number of loops of randized tests
   */

#include <recint/recint.h>

#if not defined(LOOPS)
#define LOOPS 1000
#endif

using namespace RecInt;

int main(void)
{
    rmint<STD_RECINT_SIZE> w, x, y, z;
    UDItype r;

    // Init. size = p
    RecInt::srand(limb(time(NULL)));
    ruint<STD_RECINT_SIZE> p;

    // Loop
    for (UDItype l = 1; l < LOOPS; l++) {
        do { rand(p); } while (p % 2 == 0);
        x.init_module(p);

        //------- Invert ---------

        // Not in place inv
        do { rand(x); } while (gcd(get_ruint(x), p) != 1);
        inv(y, x);
        z = x * y;
        if (z != 1) return 1;

        // In place inv
        copy(y, x);
        inv(x);
        z = x * y;
        if (z != 1) return 2;

        // Not in place inv with UDItype
        do { r = UDItype(rand()); } while (gcd(ruint<STD_RECINT_SIZE>(r), p) != 1);
        inv(x, r);
        z = x * r;
        if (z != 1 && z != 0) return 3;

        //------- Division ---------

        // Not in place div
        rand(y);
        div(z, y, x);
        z *= x;
        if (z != y && z != 0) return 4;

        // In place div
        rand(y); copy(z, y);
        div(y, x);
        y *= x;
        if (z != y && y != 0) return 5;

        // Not in place inv with UDItype
        div(z, x, r);
        z *= r;
        if (z != x && z != 0) return 6;

        // Not in place inv with UDItype
        copy(z, x);
        div(x, r);
        x *= r;
        if (z != x && x != 0) return 7;
    }

    return 0;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
