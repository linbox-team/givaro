// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

/*! @file examples/Integer/RSA_breaking.C
 * @ingroup examples
 * @ingroup integers
 * @example examples/Integer/RSA_breaking.C
 * @brief NO DOC
 */
#include <iostream>
using namespace std;
#define GIVARO_LENSTRA
#include "givaro/givintrsa.h"
#include "givaro/givtimer.h"


using namespace Givaro;




int main(int argc, char** argv)
{
    Timer tim;
    tim.clear();

    IntRSADom<>::Element m,k,u;
    if (argc > 1)
        m = IntRSADom<>::Element( argv[1] );
    else
        cin >> m;
    if (argc > 2)
        k = IntRSADom<>::Element( argv[2] );
    else
        cin >> k;

    IntRSADom<> IR(m,k);
    tim.start();
    IR.point_break(u);
    tim.stop();

    /* For a factored output :

       IR.write( cerr << "m=pq: ", IR.getm()) ;
       IR.write( cerr << ", cipher k: " , IR.getk()) << "   ----> decipering key: " ;
       IR.write( cout, u ) << endl;

*/

    // Unfactored output
    cerr << "n=pq: " << IR.getn() << ", cipher key: " << IR.gete() << "   ----> decipering key: ";
    cout << u << endl;



    cerr << tim << endl;

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
