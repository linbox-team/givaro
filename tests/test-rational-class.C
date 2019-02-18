// Copyright(c)'2018 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

#include "./rational/givrational.h"

using namespace Givaro ;

#define TEST_ASSERT(cond) \
if(!(cond)) { \
    std::cout << #cond << " failed!" << std::endl; \
    return -1; \
}

#define TEST_EQUAL(a, b) \
if(!((a) == (b))) { \
    std::cout << #a << " == " << #b << " failed! " << std::endl; \
    std::cout << (a) << " == " << (b) << std::endl; \
    return -1; \
}

int main(void)
{
    Rational zero(0);
    Rational one(1);
    Rational half(1, 2);
    Rational third(1, 3);
    Rational twoThird(2, 3);
    Rational fourSixth(4, 6);

    // ----- Basic arithmetics and auto-reduce

    TEST_ASSERT(zero != half);
    TEST_EQUAL(zero - half, -half);
    TEST_EQUAL(third + third, twoThird);
    TEST_EQUAL(2 * third, twoThird);
    TEST_EQUAL(fourSixth, twoThird);

    TEST_EQUAL(-third, Rational(-1, 3));
    TEST_EQUAL(-third, Rational(1, -3));

    // ----- round/floor/ceil

    TEST_EQUAL(round(one), 1);
    TEST_EQUAL(floor(one), 1);
    TEST_EQUAL(ceil(one), 1);

    TEST_EQUAL(round(-one), -1);
    TEST_EQUAL(floor(-one), -1);
    TEST_EQUAL(ceil(-one), -1);

    TEST_EQUAL(round(third), 0);
    TEST_EQUAL(floor(third), 0);
    TEST_EQUAL(ceil(third), 1);

    TEST_EQUAL(round(-third), 0);
    TEST_EQUAL(floor(-third), -1);
    TEST_EQUAL(ceil(-third), 0);

    TEST_EQUAL(round(twoThird), 1);
    TEST_EQUAL(floor(twoThird), 0);
    TEST_EQUAL(ceil(twoThird), 1);

    TEST_EQUAL(round(-twoThird), -1);
    TEST_EQUAL(floor(-twoThird), -1);
    TEST_EQUAL(ceil(-twoThird), 0);

    TEST_EQUAL(round(half), 1);
    TEST_EQUAL(floor(half), 0);
    TEST_EQUAL(ceil(half), 1);

    TEST_EQUAL(round(-half), -1);
    TEST_EQUAL(floor(-half), -1);
    TEST_EQUAL(ceil(-half), 0);

    // ----- Issue #74

    Rational tmp(0);
    tmp -= third;
    TEST_ASSERT(tmp.deno() >= 0);
    TEST_EQUAL(tmp, -third);

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
