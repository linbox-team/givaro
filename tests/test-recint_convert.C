/* test_rand.cpp - Randomness of RecInt library test file

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
#include <recint/ruconvert.h>

using namespace RecInt;

int main(void)
{
    ruint<STD_RECINT_SIZE> a, ca, p;
    rmint<STD_RECINT_SIZE> b, cb;
    mpz_class ga, gb;

    // Init.
    RecInt::srand(limb(time(NULL)));
    do { rand(p); } while(p % 2 == 0);
    b.init_module(p);

    // With ruint
    rand(a);
    ruint_to_mpz(ga, a);
    mpz_to_ruint(ca, ga);
    if (a != ca) return 1;

    fill_with_1(a);
    ruint_to_mpz(ga, a);
    mpz_to_ruint(ca, ga);
    if (a != ca) return 2;

    // With rmint
    rand(b);
    rmint_to_mpz(gb, b);
    mpz_to_rmint(cb, gb);
    if (b != cb) return 8;

    b = p - 1;
    rmint_to_mpz(gb, b);
    mpz_to_rmint(cb, gb);
    if (b != cb) return 9;

    return 0;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
