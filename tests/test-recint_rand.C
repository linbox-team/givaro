/* test_rand.cpp - Randomness of RecInt library test file

   Return value.
   0    No error
   != 0 Bad result for an operation

   The following constants have to be defined.
   STD_RECINT_SIZE     size of recint (> 5)
   LOOPS           number of loops of randomized tests
   */

#include <recint/recint.h>

#if not defined(LOOPS)
#define LOOPS 10000
#endif

using namespace RecInt;

int main(void)
{
    ruint<STD_RECINT_SIZE> x, p;
    rmint<STD_RECINT_SIZE> w;
    limb su = 0, tu = 0, uu = 0, vu = 0, tm = 0, lb = 0;

    // Init.
    RecInt::srand(limb(time(NULL)));
    do { rand(p); } while(p % 2 == 0);
    w.init_module(p);

    // Loop
    for (UDItype l = 1; l < LOOPS; l++) {
        // ruint
        rand(x);
        tu += limb(x & 1);
        uu += bool(x & __RECINT_MAXPOWTWO);

        lb = get_limb(x, NBLIMB<STD_RECINT_SIZE>::value-1);
        su += lb & 1;
        vu += bool(lb & __RECINT_MAXPOWTWO);

        // rmint
        rand(w);
        tm += limb(w.Value & 1);

        if (w.Value > p) return -1; // No use of module
    }

    // These tests mean that the same bit has been set to 0 a hundred times
    // (probability is < 2^100) - if so, must be an error of rand function (same part not set)
    if (LOOPS >= 100 && tu == 0) return 1;
    if (LOOPS >= 100 && tm == 0) return 2;

    if (LOOPS >= 100 && uu == 0) return 3;

    if (LOOPS >= 100 && su == 0) return 4;
    if (LOOPS >= 100 && vu == 0) return 5;

    return 0;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
