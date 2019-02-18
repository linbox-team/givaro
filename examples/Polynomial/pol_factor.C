// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
/*! @file examples/Polynomial/pol_factor.C
 * @ingroup examples
 * @ingroup polynomials
 * @example examples/Polynomial/pol_factor.C
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
    GFqDom<int64_t>::Residu_t MOD;
    if (argc > 1)
        MOD = (GFqDom<int64_t>::Residu_t) atoi(argv[1]);
    else
        std::cin >> MOD;
    uint64_t expo = 1;
    if (argc > 2) expo = (uint64_t)atoi(argv[2]);

    Timer tim2; tim2.clear(); tim2.start();
    GFqDom<int64_t> F(MOD, expo);
    tim2.stop();
    std::cerr<<"init field -> "<<tim2.usertime()<<std::endl;
    Poly1FactorDom<GFqDom<int64_t>, Dense> FD(F,Indeter("X"));
    typedef Poly1FactorDom<GFqDom<int64_t>, Dense>::Element Polys ;
    Polys P;
    FD.read( cin, P );
    std::vector<Polys> Lf;
    std::vector<uint64_t> Le;

    Timer tim; tim.clear(); tim.start();
    FD.CZfactor(Lf, Le, P);
    tim.stop();

    FD.write( cout, P ) << " is 1";
    std::vector<uint64_t>::const_iterator e = Le.begin();
    for(std::vector<Polys>::const_iterator i = Lf.begin(); i != Lf.end(); ++i, ++e) {
        FD.write(cout << " * (", *i) << ")";
        if (*e > 1) cout << "^" << *e;
    }
    std::cout << std::endl;
    std::cerr << tim << std::endl;

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
