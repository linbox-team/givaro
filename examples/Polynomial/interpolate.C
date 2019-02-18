// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
/*! @file examples/Polynomial/interpolate.C
 * @ingroup examples
 * @ingroup polynomials
 * @example examples/Polynomial/interpolate.C
 * @brief NO DOC
 */

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <givaro/modular.h>
#include <givaro/givpoly1factor.h>
#include <givaro/givtimer.h>
#include <givaro/givinterp.h>

using namespace std;

using namespace Givaro;




int main(int argc, char** argv)
{
    typedef Modular<int32_t>::Residu_t UT ;
    UT MOD;
    if (argc > 2)
        MOD = (UT) atoi(argv[2]);
    else
        std::cin >> MOD;

    Modular<int32_t> F(MOD);


    Interpolation< Modular<int32_t> > FD(F,Indeter("X"));
    Interpolation< Modular<int32_t> >::Element nouv, prec;
    int EarlyTerm = 0, Bound = 5;

    Modular<int32_t>::Element x, f;

    std::ifstream input (argv[1]);
    F.read(input, x);
    F.read(input, f);
    FD(x,f);
    prec = FD.interpolator();

    Timer tim; tim.clear(); tim.start();
    while(! F.read(input, x).eof() ) {
        F.read(input, f);
        FD(x,f);
        nouv = FD.interpolator();
        if (FD.areEqual(nouv, prec)) {
            if (++EarlyTerm > Bound) { std::cerr << "EarlyTerminated" << std::endl ; break; }
        } else
            EarlyTerm = 0;
        prec = nouv;
    }
    tim.stop();

    FD.write( std::cout, FD.interpolator()  ) << std::endl;
    std::cerr << tim << std::endl;

    return 0;
}




/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
