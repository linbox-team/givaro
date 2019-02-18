// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

/*! @file examples/Integer/ispower.C
 * @ingroup examples
 * @ingroup integers
 * @example examples/Integer/isprime.C
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


    IntPrimeDom IP;
    IntPrimeDom::Element m;
    if (argc > 1) m = Integer(argv[1]);
    unsigned int r = argc > 2 ? (unsigned int)atoi(argv[2]) : 5;

    Timer tim; tim.clear(); tim.start();
    bool a = IP.isprime(m,(int)r);
    tim.stop();
    cout << (a?"true":"false") << endl;
    cerr << tim << endl;

    //  Givaro::End();

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
