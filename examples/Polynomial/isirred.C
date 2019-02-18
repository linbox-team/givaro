// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
/*! @file examples/Polynomial/isirred.C
 * @ingroup examples
 * @ingroup polynomials
 * @example examples/Polynomial/isirred.C
 * @brief NO DOC
 */

#include <iostream>
#include <stdlib.h>
#include <givaro/gfq.h>
#include <givaro/givpoly1factor.h>
#include <givaro/givtimer.h>

using namespace std;

using namespace Givaro;




int main(int argc, char** argv)
{
    typedef GFqDom<int64_t>::Residu_t UT ;
    UT MOD;
    if (argc > 1)
        MOD =(UT) (atoi(argv[1]));
    else
        std::cin >> MOD;
    uint64_t expo = 1;
    if (argc > 2) expo = (uint64_t)atoi(argv[2]);

    GFqDom<int64_t> F(MOD, expo);

    Poly1FactorDom<GFqDom<int64_t>, Dense> FD(F,Indeter("X"));
    Poly1FactorDom<GFqDom<int64_t>, Dense>::Element P;
    FD.read( cin, P );

    Timer tim; tim.clear(); tim.start();
    bool f = FD.is_irreducible( P );
    tim.stop();

    F.write( FD.write( cout, P ) << " is " << (f?"":"not ") << "irreducible in " ) << endl;
    // std::cout << f << std::endl;
    std::cerr << tim << std::endl;

    return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
