/* rmint_neg.cpp - Negative numbers apprehension of RecInt library test file

   Return value.
   0    No error
   != 0 Bad result for an operation
   */

#include <recint/recint.h>

#if not defined(LOOPS)
#define LOOPS 10000
#endif

using namespace RecInt;

int main(void)
{
    rmint<STD_RECINT_SIZE> a, b, c;
    ruint<STD_RECINT_SIZE> p;
    DItype r;

    // Init.
    RecInt::srand(limb(time(NULL)));

    // Loop
    for (UDItype i = 1; i < LOOPS; i++) {
        do { rand(p); } while(p % 2 == 0);
        a.init_module(p);

        r = DItype(rand_gen());
        a = -r;
        b = r;
        if (a != -b) return 1;

        rand(b);
        c = p - 1; // > 0 (should equivalent to c = -1;
        a = b - c; // <= 0 (must be reset to 0..(p-1) range)
        if (a != b + 1) return 2;
    }

    return 0;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
