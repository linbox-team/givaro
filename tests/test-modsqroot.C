// =================================================================== //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <20 Jan 11 11:04:19 Jean-Guillaume.Dumas@imag.fr>
// Givaro : Modular square roots
// =================================================================== //

/*! @file examples/Integer/ModularSquareRoot.C
 * @ingroup examples
 * @ingroup integers
 * @example examples/Integer/ModularSquareRoot.C
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
    unsigned long seed = (argc>3?atoi(argv[3]):BaseTimer::seed ());
    int failures = 0;
//     std::cerr << "Seed: " << seed << std::endl;
    Integer::seeding (seed);
    Integer a,n;
    Integer ThreeToHundred;
    pow(ThreeToHundred,Integer(3),100UL);

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

    if (failures > 0) std::cerr << "test-modsqroot: " << failures << " failures." << std::endl;

    return failures;
}
