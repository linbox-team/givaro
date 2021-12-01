/* ruint_fiddling.cpp - Bits manipulation of RecInt library test file

   Return value.
   0    No error
   != 0 Bad result for an operation

   The following constants have to be defined.
   LOOPS           number of loops of randized tests
   */

#include <random>
#include <recint/ruint.h>

#if not defined(LOOPS)
#define LOOPS 10000
#endif

using namespace RecInt;

int main(void)
{
    ruint<8> x, y, z, a;
    std::mt19937_64 rand64;
    UDItype rUDI, rUDI2, rUDI3;
    limb rlimb;

    RecInt::srand(limb(time(NULL)));
    rand64.seed(limb(time(NULL)));

    // Loop
    for (UDItype l = 1; l < LOOPS; l++) {
        // Complementary and bitwise xor in place
        rand(x);
        fill_with_1(y);
        y ^= x;
        x = ~x;
        if (x != y) return 1;

        // Bitwise and in place
        rand(x);
        z = x;
        fill_with_1(y);
        x &= y;
        if (x != z) return 2;

        y &= x;
        if (y != z) return 2;

        fill_with_1(y);
        y &= 1;
        if (y != 1) return 2;

        fill_with_1(y);
        y &= __RECINT_MAXPOWTWO;
        if (y != __RECINT_MAXPOWTWO) return 2;

        fill_with_1(y);
        rUDI = rand64();
        y &= rUDI;
        if (y != rUDI) return 2;

        // Bitwise or in place
        rand(x);
        fill_with_1(y);
        z = y;
        x |= y;
        if (x != z) return 3;

        rand(x);
        y |= x;
        if (y != z) return 3;

        reset(y);
        y |= 1;
        if (y != 1) return 3;

        reset(y);
        y |= __RECINT_MAXPOWTWO;
        if (y != __RECINT_MAXPOWTWO) return 3;

        // Bitwise and
        rand(x);
        rUDI = rand64();
        rUDI2 = x & rUDI;
        rlimb = get_limb(x, 0);
        rUDI3 = rlimb;
        rUDI3 &= rUDI;
        if (rUDI2 != rUDI3) return 4;

        // Bitwise or
        rand(x);
        rUDI = rand64();
        z = x | rUDI;
        rlimb = get_limb(x, 0);
        rUDI3 = rlimb;
        rUDI3 |= rUDI;
        set_limb(x, limb(rUDI3), 0);
        if (x != z) return 5;

        // Bitwise xor
        rand(x);
        rUDI = rand64();
        z = x ^ rUDI;
        rlimb = get_limb(x, 0);
        rUDI3 = rlimb;
        rUDI3 ^= rUDI;
        set_limb(x, limb(rUDI3), 0);
        if (x != z) return 6;

        // Classic test
        rand(x); rand(y);
        z = x^y;
        a = ((~x) & y) | (x & (~y));
        if (a != z) return 7;

        rlimb = rand64() & 31;
        reset(x);
        set_highest_word(x, __RECINT_MAXPOWTWO>>rlimb );
        x <<= rlimb;
        max_pow_two(y);
        if (x != y) {
            std::cerr << x << " != " << y << std::endl;
            return 10;
        }
    }

    // Extra functions
    if (highest_bit(max_pow_two(x)) != true) return 8;
    if (lowest_bit(~x) != true) return 8;

    set_highest_word(x,__RECINT_MAXPOWTWO);
    max_pow_two(y);
    if (x != y) return 9;
    

    return 0;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
