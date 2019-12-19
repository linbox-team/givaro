// Copyright(c)'1994-2025 by The Givaro group
// This file is part of Givaro.
// BB: adapted,enhanced from examples/FiniteFields/ff-arith.C
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

#include <iostream>
#include "test-fieldarith.h"
#ifdef __GIVARO_DEBUG
#define GIVARO_RATRECON_DEBUG
#endif

#include <givaro/qfield.h>

template<>
bool invertible(const QField<Rational>& Q, const Rational& a)
{
    return !Q.isZero(a);
}

int main(int argc, char ** argv)
{
    int seed = int (argc>1?atoi(argv[1]):BaseTimer::seed());
#ifdef __GIVARO_DEBUG
    std::cerr << "seed: " << seed << std::endl;
#endif
    Integer::seeding((uint64_t)seed);
        // Rational numbers

    TEST_SPECIFIC(QField<Rational>, QQ, 0);
    return 0;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
