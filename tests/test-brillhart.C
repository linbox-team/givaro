// =================================================================== //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <29 Jun 19 08:40:29 Jean-Guillaume.Dumas@imag.fr>
// Givaro : Modular square roots
// =================================================================== //

/*! @file tests/test-modsqroot.C
 * @ingroup examples
 * @ingroup integers
 * @see examples/Integer/ModularSquareRoot.C
 * @brief NO DOC
 */
#include <iostream>
#include <stdlib.h>
#include <givaro/givintsqrootmod.h>
#include <givaro/givtimer.h>

using namespace Givaro;
IntSqrtModDom<> ISM;

bool TestBrillart(const Integer& p) {
    Integer a,b;
    ISM.Brillhart(a,b,p);
    return ISM.areEqual(a*a+b*b,p);
}

bool TestSoSmP(const Integer& k, const Integer& p) {
    Integer a,b;
    ISM.sumofsquaresmodprime(a,b,k,p);
    return ISM.isZero( (a*a+b*b-k) % p );
}



int main(int argc, char** argv) {
    int nbtests = (argc>1?atoi(argv[1]):1000);
    int sizes = (argc>2?atoi(argv[2]):100);
    unsigned long seed = (unsigned long)(argc>3?(unsigned long)atoi(argv[3]):(unsigned long)BaseTimer::seed ());
#ifdef GIVARO_DEBUG
    std::cerr << "seed: " << seed << std::endl;
#endif
    int failures = 0;


    Integer::seeding (seed);
    Integer a,n;

    if (! TestBrillart(Integer(10006721))) ++failures;

    for(int i=0; i<nbtests; ++i) {

        Integer::random(n,sizes);
        do {
            ISM.nextprimein(n);
        } while( ISM.mod(a,n,4U) != 1);

        if (! TestBrillart(n)) ++failures;
    }

    if (failures > 0) std::cerr << "Brillhart: " << failures << " failures." << std::endl;

    Integer k(ISM.mOne),b; n=23;
    ISM.sumofsquaresmodprime(a,b,k,n);
    if (!ISM.isZero( (a*a+b*b-k) % n )) ++failures;
    else
        std::clog << a << '*' << a << '+' << b << '*' << b << '=' << k << '%' << n << std::endl;

    for(int i=0; i<nbtests; ++i) {

        Integer::random(a,sizes);
        Integer::random(n,sizes);
        ISM.nextprimein(n);
        ISM.modin(a,n);
        if (! TestSoSmP(a,n)) ++failures;
    }

    if (failures > 0) std::cerr << "Modular SoS: " << failures << " failures." << std::endl;

    return failures;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
