// =================================================================== //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <18 Jan 11 14:37:04 Jean-Guillaume.Dumas@imag.fr>
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

int main(int argc, char** argv) {
    
    int nbtests = (argc>1?atoi(argv[1]):100);
    int sizes = (argc>2?atoi(argv[2]):10);
    int failures = 0;
    Integer::seeding (BaseTimer::seed ());
    IntSqrtModDom<> ISM;
    Integer r,x,a,b,n(1UL);

    for(int i=0; i<nbtests; ++i) {
        
        Integer::random(a,sizes);
        while (n<=1) Integer::nonzerorandom(n,sizes);

        ISM.modin(a,n);
        
        ISM.modin(ISM.mul(b,a,a),n); // b is now a square mod n
        ISM.sqrootmod(x,b,n);
        if (! ISM.areEqual( ISM.modin(ISM.mul(r,x,x),n), b) ) {
            
            std::cerr << "n: " << n << std::endl;
            std::cerr << "a: " << a << std::endl;
            std::cerr << "x: " << x << std::endl;
            std::cerr << "r: " << r << std::endl;
            
            ++failures;
        }
        
    }
     
    return failures;
}
