// =================================================================== //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <28 Jun 19 18:15:50 Jean-Guillaume.Dumas@imag.fr>
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

bool TestSrqtMod(const Integer& z, const Integer& n) {
    Integer a(z),b,x,r;
    ISM.modin(a,n);
    ISM.modin(ISM.mul(b,a,a),n); // b is now a square mod n
    ISM.sqrootmod(x,b,n);
    if (ISM.areEqual( ISM.modin(ISM.mul(r,x,x),n), b) ) {
        return true;
    } else {
        std::cerr << "TestSrqtMod ERROR" << std::endl;
        std::cerr << "n:= " << n << ';' << std::endl;
        std::cerr << "a:= " << a << ';' << std::endl;
        std::cerr << "x:= " << x << ';' << std::endl;
        std::cerr << "r:= " << r << ';' << std::endl;

        return false;
    }
}

int main(int argc, char** argv) {
    int nbtests = (argc>1?atoi(argv[1]):100);
    int sizes = (argc>2?atoi(argv[2]):10);
    unsigned long seed = (unsigned long)(argc>3?(unsigned long)atoi(argv[3]):(unsigned long)BaseTimer::seed ());
#ifdef __GIVARO_DEBUG
    std::cerr << "seed: " << seed << std::endl;
#endif
    int failures = 0;


    //     std::cerr << "Seed: " << seed << std::endl;
    Integer::seeding (seed);
    Integer a,n;
    Integer ThreeToHundred;
    pow(ThreeToHundred,Integer(3),100U);

    for(int i=0; i<nbtests; ++i) {

        Integer::random(a,sizes);
        Integer::nonzerorandom(n,sizes);
        while(n<=1) Integer::nonzerorandom(n,sizes);

        if (! TestSrqtMod(a,n)) ++failures;

        n <<= (129);
        n *= ThreeToHundred;

        //         ISM.write(std::cerr, n) << std::endl;

        if (! TestSrqtMod(a,n)) ++failures;


    }

    if (failures > 0) std::cerr << "modsqroot: " << failures << " failures." << std::endl;

    return failures;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
