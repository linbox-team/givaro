// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

/*! @file examples/Integer/ifactor_lenstra.C
 * @ingroup examples
 * @ingroup integers
 * @example examples/Integer/ifactor_lenstra.C
 * @brief NO DOC
 */
#define GIVARO_LENSTRA
#include <iostream>
using namespace std;
#include <givaro/givinit.h>


#include <givaro/givintfactor.h>
#include <givaro/givtimer.h>

using namespace Givaro;



int main(int argc, char** argv)
{
    IntFactorDom<> IP;
#ifndef __GIVARO_GMP_NO_CXX
    IP.seeding();
    // std::cerr << "Seeding..." << std::endl;
#endif
    Integer m;
    if (argc > 1)
        m = Integer(argv[1]);
    else
        cin >> m;
    if (IP.islt(m,0) ) {
        cerr << "-";
        IP.negin(m);
    }
    if (IP.islt(m,4))
        IP.write(cerr,m) << endl;
    else {
        Timer tim; tim.clear(); tim.start();
        IP.write(cerr,m) << endl;
        tim.stop();
        cerr << tim << endl;
    }
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
