// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

/*! @file examples/Integer/isproot.C
 * @ingroup examples
 * @ingroup integers
 * @example examples/Integer/isproot.C
 * @brief NO DOC
 */
#include <iostream>
#include <givaro/givintnumtheo.h>
#include <givaro/givtimer.h>

using namespace Givaro;



//----------------------------------------------------------
// Deterministic, non-polynomial (factor b-1 for the order),
// test for primitive roots.
//----------------------------------------------------------

int main(int argc, char** argv)
{
    IntNumTheoDom<> IP;
    IntNumTheoDom<>::Element a,b;
    if (argc > 1) a = Integer(argv[1]); else std::cin >> a;
    if (argc > 2) b = Integer(argv[2]); else std::cin >> b;

    Timer tim; tim.clear(); tim.start();
    bool f = IP.is_prim_root(a,b);
    tim.stop();
    std::cout << f << std::endl;
    std::cerr << tim << std::endl;

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
