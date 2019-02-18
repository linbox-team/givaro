/* rmint_constructors.cpp - class rmint constructors of RecInt test file

   Return value.
   0    No error
   != 0 Bad result for an operation
   */

#include <recint/recint.h>

using namespace RecInt;

int main(void)
{
    ruint<10> p10;
    ruint<__RECINT_LIMB_SIZE> pLS;

    RecInt::srand(limb(time(NULL)));
    rand(p10); if ((p10 % 2) == 0) p10++;
    rand(pLS); if ((pLS % 2) == 0) pLS++;

    rmint<10>::init_module(p10); // Initializing module for rmint<10>
    rmint<__RECINT_LIMB_SIZE>::init_module(pLS); // The same but for rmint<__RECINT_LIMB_SIZE>

    // Large ruint, small value
    rmint<10> x1(1);
    rmint<10> x2(UDItype(1));
    rmint<10> x3(USItype(1));
    rmint<10> x4(DItype(1));
    rmint<10> x5(SItype(1));

    if (x1 != x2) return 1;
    if (x1 != x3) return 2;
    if (x1 != x4) return 3;
    if (x1 != x5) return 4;
    if (x1 != 1)  return 5;

    // Small ruint, small value
    rmint<__RECINT_LIMB_SIZE> y1(1);
    rmint<__RECINT_LIMB_SIZE> y2(UDItype(1));
    rmint<__RECINT_LIMB_SIZE> y3(USItype(1));
    rmint<__RECINT_LIMB_SIZE> y4(DItype(1));
    rmint<__RECINT_LIMB_SIZE> y5(SItype(1));

    if (y1 != y2) return 6;
    if (y1 != y3) return 7;
    if (y1 != y4) return 8;
    if (y1 != y5) return 9;
    if (y1 != 1)  return 10;

    // Large ruint, negative value
    rmint<10> v1(DItype(-1));
    rmint<10> v2(SItype(-1));
    rmint<10> v3(p10 - 1);

    if (v1 != v2) return 11;
    if (v1 != v3) return 12;

    // Small ruint, negative value
    rmint<__RECINT_LIMB_SIZE> w1(DItype(-1));
    rmint<__RECINT_LIMB_SIZE> w2(SItype(-1));
    rmint<__RECINT_LIMB_SIZE> w3(pLS - 1);

    if (w1 != w2) return 13;
    if (w1 != w3) return 14;

    // Copy
    rmint<10> xi1(x1), vi1(v1);
    rmint<__RECINT_LIMB_SIZE> yi1(y1), wi1(w1);

    if (xi1 != x1) return 15;
    if (vi1 != v1) return 16;
    if (yi1 != y1) return 17;
    if (wi1 != w1) return 18;

    // rmint from ruint
    // check that it is correctly reduced
    rmint<10> pi(p10);
    rmint<__RECINT_LIMB_SIZE> li(pLS);

    if (pi != 0) return 19;
    if (li != 0) return 20;

    return 0;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
