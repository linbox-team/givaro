// =================================================================== //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <15 Jul 19 10:30:33 Jean-Guillaume.Dumas@imag.fr>
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

Givaro::Timer sosb, sosl;

bool TestSoSmP(const Integer& k, const Integer& p) {
    bool pass(true);

    Integer a,b;
    Givaro::Timer chrono; chrono.start();
    ISM.sumofsquaresmodprimeNoERH(a,b,k,p);
    chrono.stop();
    sosb += chrono;
    pass &= ISM.isZero( (a*a+b*b-k) % p );

    chrono.start();
    ISM.sumofsquaresmodprime(a,b,k,p);
    chrono.stop();
    sosl += chrono;
    pass &= ISM.isZero( (a*a+b*b-k) % p );

    return pass;
}

int main(int argc, char** argv) {
    uint64_t nbtests = (argc>1?atoi(argv[1]):1000);
    uint64_t sizes = (argc>2?atoi(argv[2]):100);
    uint64_t seed = (uint64_t)(argc>3?(uint64_t)atoi(argv[3]):(uint64_t)BaseTimer::seed ());
#ifdef GIVARO_DEBUG
    std::cerr << "seed: " << seed << std::endl;
#endif
    uint64_t failures = 0;

    sosb.clear();
    sosl.clear();

    Integer::seeding (seed);
    Integer a,n;

    if (! TestBrillart(Integer(10006721))) ++failures;

    for(uint64_t i=0; i<nbtests; ++i) {

        Integer::random(n,sizes);
        do {
            ISM.nextprimein(n);
        } while( ISM.mod(a,n,4U) != 1);

        if (! TestBrillart(n)) ++failures;
    }

    if (failures > 0) std::cerr << "Brillhart: " << failures << " failures." << std::endl;

    Integer k(ISM.mOne),b; n=23;
    if (! TestSoSmP(k,n) ) ++failures;

    for(uint64_t i=0; i<nbtests; ++i) {

        Integer::random(a,sizes);
        Integer::random(n,sizes);
        ISM.nextprimein(n);
        ISM.modin(a,n);
        if (! TestSoSmP(a,n)) ++failures;
    }

    if (failures > 0) std::cerr << "Modular SoS: " << failures << " failures." << std::endl;


    std::clog << "SOSB: " << sosb << std::endl;
    std::clog << "SOSL: " << sosl << std::endl;


    return failures;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
