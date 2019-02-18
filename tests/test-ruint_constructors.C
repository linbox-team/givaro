/* test_operators.cpp - Arithmetic constructors of RecInt test file

   Return value.
   0    No error
   != 0 Bad result for an operation
   */

#include <recint/ruint.h>

using namespace RecInt;

template<size_t SIZEOFTEST=STD_RECINT_SIZE>
int mainTest(void)
{
    // Large ruint, small value
    ruint<SIZEOFTEST> x1(bool(true));
    ruint<SIZEOFTEST> x2(UDItype(1));
    ruint<SIZEOFTEST> x3(USItype(1));
    ruint<SIZEOFTEST> x4(DItype(1));
    ruint<SIZEOFTEST> x5(SItype(1));

    if (x1 != x2) return 1;
    if (x1 != x3) return 2;
    if (x1 != x4) return 3;
    if (x1 != x5) return 4;

    // Small ruint, small value
    ruint<__RECINT_LIMB_SIZE> y1(bool(true));
    ruint<__RECINT_LIMB_SIZE> y2(UDItype(1));
    ruint<__RECINT_LIMB_SIZE> y3(USItype(1));
    ruint<__RECINT_LIMB_SIZE> y4(DItype(1));
    ruint<__RECINT_LIMB_SIZE> y5(SItype(1));

    if (y1 != y2) return 5;
    if (y1 != y3) return 6;
    if (y1 != y4) return 7;
    if (y1 != y5) return 8;

    // Large ruint, negative value
    ruint<SIZEOFTEST> v1(DItype(-1));
    ruint<SIZEOFTEST> v2(SItype(-1));
    ruint<SIZEOFTEST> v3;

    fill_with_1(v3);
    if (v1 != v2) return 9;
    if (v1 != v3) return 10;

    // Small ruint, negative value
    ruint<__RECINT_LIMB_SIZE> w1(DItype(-1));
    ruint<__RECINT_LIMB_SIZE> w2(SItype(-1));

    if (w1 != w2) return 11;
    if (w1 != __RECINT_MINUSONE) return 12;

    // Copy test
    ruint<SIZEOFTEST> xi1(x1), vi1(v1);
    ruint<__RECINT_LIMB_SIZE> yi1(y1), wi1(w1);

    if (xi1 != x1) return 13;
    if (vi1 != v1) return 14;
    if (yi1 != y1) return 15;
    if (wi1 != w1) return 16;

    return 0;
}


int main(void) {
    int rtm = mainTest<>();
    rtm += mainTest<10>();

    return rtm;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
