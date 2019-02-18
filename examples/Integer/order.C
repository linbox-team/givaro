// ========================================================== //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <04 Sep 02 18:10:22 Jean-Guillaume.Dumas@imag.fr>
// ========================================================== //
/*! @file examples/Integer/order.C
 * @ingroup examples
 * @ingroup integers
 * @example examples/Integer/order.C
 * @brief NO DOC
 */
#include <iostream>
using namespace std;
#include <stdlib.h>
#include <givaro/givintnumtheo.h>
#include <givaro/givtimer.h>

using namespace Givaro;




int main(int argc, char ** argv)
{
    IntNumTheoDom<> IP;
    IntNumTheoDom<>::Element a,q,o;
    if (argc > 1) a = Integer(argv[1]); else cin >> a;
    if (argc > 2) q = Integer(argv[2]); else cin >> q;

    Timer tim; tim.clear(); tim.start();
    // Ordre de a dans GF(q)
    IP.order(o, a, q);
    tim.stop();

    cout << o << endl;
    cerr << tim << endl;

    return 0;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
