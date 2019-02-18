/* test_operators.cpp - Arithmetic operators of RecInt generic test file

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
#define LOOPS 10000
#endif

using namespace RecInt;

int main(void)
{
    ruint<STD_RECINT_SIZE> x, y, z, s;
    mpz_class size, gx, gy, gz, gs, gcmp;
    USItype r;

    // Init. size = 2 ^ (2 ^ STD_RECINT_SIZE)
    mpz_ui_pow_ui(size.get_mpz_t(), 2, STD_RECINT_SIZE);
    mpz_ui_pow_ui(size.get_mpz_t(), 2, size.get_ui());
    RecInt::srand(limb(time(NULL)));

    // Loop
    for (UDItype l = 1; l < LOOPS; l++) {
        // RecInt rand
        rand(x);
        rand(y);
        rand(z);
        ruint_to_mpz(gx, x);
        ruint_to_mpz(gy, y);
        ruint_to_mpz(gz, z);

        // Comparison
        if (x < x) return -2;
        if (x > x) return -2;
        if (x != x) return -2;

        if (x < y && gx >= gy) return -1;
        if (x > y && gx <= gy) return -1;
        if (x <= y && gx > gy) return -1;
        if (x >= y && gx < gy) return -1;
        if (x == y && gx != gy) return -1;
        if (x != y && gx == gy) return -1;

        // Increment and decrement
        y++; gy++; ruint_to_mpz(gcmp, y);
        if (gcmp != gy) return 1;
        ++y; ++gy; ruint_to_mpz(gcmp, y);
        if (gcmp != gy) return 1;

        if (x < 2) { x = 2; gx = 2; }
        x--; gx--; ruint_to_mpz(gcmp, x);
        if (gcmp != gx) return 2;
        --x; --gx; ruint_to_mpz(gcmp, x);
        if (gcmp != gx) return 2;

        // add and sub in place
        x += y; gx += gy; gx %= size; ruint_to_mpz(gcmp, x);
        if (gcmp != gx) return 3;
        r = USItype(rand()); x += r; gx += r; gx %= size; ruint_to_mpz(gcmp, x);
        if (gcmp != gx) return 3;
        rand(s); ruint_to_mpz(gs, s); s += s; gs += gs; gs %= size; ruint_to_mpz(gcmp, s);
        if (gs != gcmp) return 3;

        if (z > x) {
            z -= x; gz -= gx; ruint_to_mpz(gcmp, z);
            if (gcmp != gz) return 4;
        } else {
            x -= z; gx -= gz; ruint_to_mpz(gcmp, x);
            if (gcmp != gx) return 4;
        }
        while (y < r) { rand(y); r = USItype(rand()); }
        ruint_to_mpz(gy, y); y -= r; gy -= r; ruint_to_mpz(gcmp, y);
        if (gcmp != gy) return 4;
        rand(s); ruint_to_mpz(gs, s); s -= s; gs -= gs; ruint_to_mpz(gcmp, s);
        if (gs != gcmp) return 4;

        // mul and div in place
        if (z < 2) { z = 2; gz = 2; }
        x %= z; gx %= gz;
        ruint_to_mpz(gcmp, x);
        if (gcmp != gx) return 5;
        while ((r = USItype(rand())) < 2);
        z %= r; gz %= r; ruint_to_mpz(gcmp, z);
        if (gcmp != gz) return 5;
        while (s == 0) rand(s);
        ruint_to_mpz(gs, s); s %= s; gs %= gs; ruint_to_mpz(gcmp, s);
        if (gs != gcmp) return 5;

        x *= y; gx *= gy; gx %= size; ruint_to_mpz(gcmp, x);
        if (gcmp != gx) return 6;
        r = USItype(rand()); x *= r; gx *= r; gx %= size; ruint_to_mpz(gcmp, x);
        if (gcmp != gx) return 6;
        rand(s); ruint_to_mpz(gs, s); s *= s; gs *= gs; gs %= size; ruint_to_mpz(gcmp, s);
        if (gs != gcmp) return 6;

        if (y == 0) { y++; gy++; }
        z /= y; gz /= gy;
        ruint_to_mpz(gcmp, z);
        if (gcmp != gz) return 7;
        while ((r = USItype(rand())) < 2);
        y /= r; gy /= r; ruint_to_mpz(gcmp, y);
        if (gcmp != gy) return 7;
        while (s == 0) rand(s);
        ruint_to_mpz(gs, s); s /= s; gs /= gs; ruint_to_mpz(gcmp, s);
        if (gs != gcmp) return 7;

        // Refresh
        rand(x); ruint_to_mpz(gx, x);
        rand(y); ruint_to_mpz(gy, y);
        rand(z); ruint_to_mpz(gz, z);

        // + symbol
        x = y + z; gx = gy + gz; gx %= size; ruint_to_mpz(gcmp, x);
        if (gcmp != gx) return 8;
        r = USItype(rand()); x = (UDItype)r + z; gx = r + gz; gx %= size; ruint_to_mpz(gcmp, x);
        if (gcmp != gx) return 8;
        r = USItype(rand()); x = y + (UDItype)r; gx = gy + r; gx %= size; ruint_to_mpz(gcmp, x);
        if (gcmp != gx) return 8;

        // - symbol
        while (y < z) { rand(y); rand(z); }
        ruint_to_mpz(gy, y); ruint_to_mpz(gz, z); x = y - z; gx = gy - gz; ruint_to_mpz(gcmp, x);
        if (gcmp != gx) return 9;
        while (r < z) { r = USItype(rand()); z = rand(); }
        ruint_to_mpz(gz, z); x = (UDItype)r - z; gx = r - gz; ruint_to_mpz(gcmp, x);
        if (gcmp != gx) return 9;
        while (y < r) { r = USItype(rand()); rand(y); }
        ruint_to_mpz(gy, y); x = y - (UDItype)r; gx = gy - r; ruint_to_mpz(gcmp, x);
        if (gcmp != gx) return 9;

        // Refresh
        rand(x); ruint_to_mpz(gx, x);
        rand(y); ruint_to_mpz(gy, y);
        rand(z); ruint_to_mpz(gz, z);

        // % symbol
        z = x % y;
        gz = gx % gy; ruint_to_mpz(gcmp, z);
        if (gcmp != gz) return 10;
        while ((r = USItype(rand())) < 2);
        z = (x % (UDItype)r); gz = gx % r; ruint_to_mpz(gcmp, z);
        if (gcmp != gz) return 10;
        if (x % x != 0) return 10;

        // * symbol
        x = y * z; gx = gy * gz; gx %= size; ruint_to_mpz(gcmp, x);
        if (gcmp != gx) return 11;
        r = USItype(rand()); x = (UDItype)r * z; gx = r * gz; gx %= size; ruint_to_mpz(gcmp, x);
        if (gcmp != gx) return 11;
        r = USItype(rand()); x = y * (UDItype)r; gx = gy * r; gx %= size; ruint_to_mpz(gcmp, x);
        if (gcmp != gx) return 11;

        // / symbol
        z = x / y; gz = gx / gy; ruint_to_mpz(gcmp, z);
        if (gcmp != gz) return 12;
        while ((r = USItype(rand())) < 2);
        z = x / (UDItype)r; gz = gx / r; ruint_to_mpz(gcmp, z);
        if (gcmp != gz) return 12;
        if (x / x != 1) return 12;
    }

    return 0;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
