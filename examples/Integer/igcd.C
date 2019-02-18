// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

/*! @file examples/Integer/igcd.C
 * @ingroup examples
 * @ingroup integers
 * @example examples/Integer/igcd.C
 * @brief NO DOC
 */
#include <iostream>
using namespace std;
#include <stdlib.h>
#include <givaro/givintprime.h>
#include <givaro/givtimer.h>
#include <givaro/givinit.h>         // Givaro initialization


using namespace Givaro;




int main(int argc, char** argv)
{
    //  Givaro::Init(&argc, &argv);


    IntegerDom IP;
    IntegerDom::Element GG, g,a,b;
    int offset = 0;
    if (argc > ++offset) a = Integer(argv[offset]); else cin >> a;
    if (argc > ++offset) b = Integer(argv[offset]); else cin >> b;

    Timer tim; tim.clear(); tim.start();
    IP.gcd(GG,a,b);
    for ( ; argc > ++offset; ) {
        a = Integer(argv[offset]);
        IP.gcd(g, GG, a);
        GG = g;
    }
    tim.stop();
    cout << GG << endl;
    cerr << tim << endl;

    //  Givaro::End();

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
