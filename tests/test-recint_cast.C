/* recint_cast.cpp - Cast functions of RecInt library test file

   Return value.
   0    No error
   != 0 Bad result for an operation

   The following constants have to be defined.
   LOOPS           number of loops of randomized tests
   */

#define MG_DEFAULT MG_INACTIVE
#include <recint/recint.h>

#if not defined(LOOPS)
#define LOOPS 10000
#endif

using namespace RecInt;

int main(void)
{
    ruint<STD_RECINT_SIZE> a, zero(0), p;
    rmint<STD_RECINT_SIZE> m, zerom;
    rmint<STD_RECINT_SIZE, MG_ACTIVE> mg, mg2, zeromg;

    UDItype l;
    limb lb;

    // Init.
    RecInt::srand(limb(time(NULL)));

    // Loop
    for (UDItype i = 1; i < LOOPS; i++) {
        do { rand(p); } while(p % 2 == 0);
        m.init_module(p); // For basic
        mg.init_module(p); // For montgomery

        //--------- ruint -----------

        rand(a);
        lb = get_limb(a, 0);
        l = UDItype(lb);

        // bool
        if (a != 0 && !bool(a)) return 1;
        if ( (a<<33) != 0 && !bool(a<<33)) return 133;
        if ( (STD_RECINT_SIZE>6) && (a<<65) != 0 && !bool(a<<65)) return 165;
        if ( (STD_RECINT_SIZE>7) && (a<<129) != 0 && !bool(a<<129)) return 1129;
        if (bool(zero)) return 1;

        // UDItype
        if (UDItype(a) != l) return 2;
        if (USItype(a) != USItype(l)) return 3;
        //        if (DItype(a) != DItype(l)) return 4; This has changed for GMP compatibility
        //        if (SItype(a) != SItype(l)) return 5; This has changed for GMP compatibility


        //--------- rmint no Montgomery -----------

        rand(m);
        lb = get_limb(m.Value, 0);
        l = UDItype(lb);

        // bool
        if (m != 0 && !bool(m)) return 6;
        if (bool(zerom)) return 6;

        // UDItype
        if (UDItype(m) != l) return 7;
        if (USItype(m) != USItype(l)) return 8;
        //        if (DItype(m) != DItype(l)) return 9;  This has changed for GMP compatibility
        //        if (SItype(m) != SItype(l)) return 10; This has changed for GMP compatibility


        //--------- rmint with Montgomery -----------

        rand(mg);
        reduction(mg2, mg);
        lb = get_limb(mg2.Value, 0);
        l = UDItype(lb);

        // bool
        if (mg != 0 && !bool(mg)) return 11;
        if (bool(zeromg)) return 11;

        // UDItype
        if (UDItype(mg) != l) return 12;
        if (USItype(mg) != USItype(l)) return 13;
        //        if (DItype(mg) != DItype(l)) return 14; This has changed for GMP compatibility
        //        if (SItype(mg) != SItype(l)) return 15; This has changed for GMP compatibility
    }

    return 0;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
