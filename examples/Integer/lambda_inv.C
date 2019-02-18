// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

/*! @file examples/Integer/lambda_inv.C
 * @ingroup examples
 * @ingroup integers
 * @example examples/Integer/lambda_inv.C
 * @brief NO DOC
 */
#include <iostream>
#include <givaro/givintnumtheo.h>
#include <givaro/givtimer.h>

using namespace Givaro;



// Lambda function : order of a primitive invertible
//                (invertible Element of maximal orbit size)
//


int main(int argc, char** argv)
{
    IntNumTheoDom<> IP;
    IntNumTheoDom<>::Element a,pr;
    if (argc > 1) a = IntNumTheoDom<>::Element(argv[1]); else std::cin >> a;

    Timer tim; tim.clear(); tim.start();
    IP.lambda_inv(pr, a);
    tim.stop();
    IntegerDom().write( std::cout, pr ) << std::endl;
    std::cerr << tim << std::endl;

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
