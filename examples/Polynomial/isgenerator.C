// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
/*! @file examples/Polynomial/isprimitive.C
 * @ingroup examples
 * @ingroup polynomials
 * @example examples/Polynomial/isprimitive.C
 * @brief NO DOC
 */

#include <iostream>
#include <stdlib.h>
#include <givaro/gfq.h>
#include <givaro/givpoly1factor.h>
#include <givaro/givtimer.h>

// using namespace std;

using namespace Givaro;


int main(int argc, char** argv)
{
    typedef GFqDom<int64_t>::Residu_t UT;
    UT MOD;
    if (argc > 1)
        MOD = (UT)atoi(argv[1]);
    else
        std::cin >> MOD;
    uint64_t expo = 1;
    if (argc > 2) expo = (uint64_t)atoi(argv[2]);
    
    GFqDom<int64_t> F(MOD, expo);
    
    Poly1FactorDom<GFqDom<int64_t>, Dense> FD(F,Indeter("X"));
    Poly1FactorDom<GFqDom<int64_t>, Dense>::Element P, G;
    FD.read( std::cin, P );
    FD.read( std::cin, G );
    
    Timer tim; tim.clear(); tim.start();
    bool f = FD.is_prim_root(G, P );
    tim.stop();
    
    FD.write(F.write( FD.write( std::cout, G ) << " is " << (f?"":"not ") << " a generator in " ) << " defined with ", P) << std::endl;
    
        // std::cout << f << std::endl;
    std::cerr << tim << std::endl;
    
    return 0;
}
