// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
/*! @file examples/Polynomial/pol_eval.C
 * @ingroup examples
 * @ingroup polynomials
 * @example examples/Polynomial/pol_eval.C
 * @brief NO DOC
 */

#include <iostream>
#include <stdlib.h>
#include <givaro/gfq.h>
#include <givaro/givpoly1factor.h>
#include <givaro/givtimer.h>

using namespace Givaro;


using namespace std;


int main(int argc, char** argv)
{
    GFqDom<int64_t>::Residu_t MOD;
    if (argc > 1)
        MOD = (GFqDom<int64_t>::Residu_t) atoi(argv[1]);
    else
        std::cin >> MOD;
    uint64_t expo = 1;
    if (argc > 2) expo = (uint64_t)atoi(argv[2]);

    GFqDom<int64_t> F(MOD, expo);

    Poly1FactorDom<GFqDom<int64_t>, Dense> FD(F,Indeter("X"));
    Poly1FactorDom<GFqDom<int64_t>, Dense>::Element P;
    FD.read( cin, P );
    GFqDom<int64_t>::Element res, val;
    F.read( cin, val );

    Timer tim; tim.clear(); tim.start();
    FD.eval(res, P, val );
    tim.stop();

    F.write( F.write( FD.write( cout, P ) << " is ", res ) << " at ", val) << endl;
    std::cerr << tim << std::endl;

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
