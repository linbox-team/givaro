// =================================================================== //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <19 Nov 19 13:07:56 Jean-Guillaume.Dumas@imag.fr>
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

Givaro::Timer sosb, sosr, sosl, brill;

bool TestBrillart(const Integer& p) {
    Integer a,b;
    Givaro::Timer chrono; chrono.start();
    ISM.Brillhart(a,b,p);
    chrono.stop();
    brill += chrono;
    return ISM.areEqual(a*a+b*b,p);
}

bool TestSoSmP(const Integer& k, const Integer& p) {
    bool pass(true);

    Integer a,b;
    Givaro::Timer chrono; chrono.start();
    ISM.sumofsquaresmodprimeNoERH(a,b,k,p);
    chrono.stop();
    sosb += chrono;
    pass &= ISM.isZero( (a*a+b*b-k) % p );

    chrono.start();
    ISM.sumofsquaresmodprimeDeterministic(a,b,k,p);
    chrono.stop();
    sosl += chrono;
    pass &= ISM.isZero( (a*a+b*b-k) % p );

    chrono.start();
    ISM.sumofsquaresmodprimeMonteCarlo(a,b,k,p);
    chrono.stop();
    sosr += chrono;
    pass &= ISM.isZero( (a*a+b*b-k) % p );

    return pass;
}

int main(int argc, char** argv) {
    uint64_t nbtests = (argc>1?atoi(argv[1]):1000);
    uint64_t sizes = (argc>2?atoi(argv[2]):100);
    uint64_t seed = (uint64_t)(argc>3?(uint64_t)atoi(argv[3]):(uint64_t)BaseTimer::seed ());
#ifdef __GIVARO_DEBUG
    std::cerr << "seed: " << seed << std::endl;
#endif
    uint64_t failures = 0;

    sosb.clear();
    sosr.clear();
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
    std::clog << "BRILL: " << brill << std::endl;

    Integer k(ISM.mOne),b; n=23;
    if (! TestSoSmP(k,n) ) ++failures;

    uint64_t nbp( uint64_t( Givaro::sqrt(nbtests) ) );
    if (nbp<0) nbp = 1;
    uint64_t nba( ceil(nbtests/nbp) );

    for(uint64_t i=0; i<nbp; ++i) {
        Integer::random(n,sizes);
        ISM.nextprimein(n);
        for(uint64_t j=0; j<nba; ++j) {
            Integer::random(a,sizes);
            ISM.modin(a,n);
            if (! TestSoSmP(a,n)) ++failures;
        }
    }

    if (failures > 0) std::cerr << "Modular SoS: " << failures << " failures." << std::endl;
    else {
        std::clog << "[SOSB] PASSED: " << sosb << std::endl;
        std::clog << "[SOSL] PASSED: " << sosl << std::endl;
        std::clog << "[SOSR] PASSED: " << sosr << std::endl;
    }

    return failures;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
