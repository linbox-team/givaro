// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

/*! @file examples/Integer/RSA_encipher.C
 * @ingroup examples
 * @ingroup integers
 * @example examples/Integer/RSA_encipher.C
 * @brief NO DOC
 */
#include <iostream>
#include <fstream>
#include "givaro/givintrsa.h"
#include "givaro/givrandom.h"
#include "givaro/givtimer.h"

// RSA, in CBC mode, enciphering of files

using namespace Givaro;




int main(int argc, char** argv)
{
    Timer tim;
    tim.clear();

    IntRSADom<GivRandom>::Rep n,e;
    if (argc > 3)
        n = Integer( argv[3] );
    else
        std::cin >> n;
    if (argc > 4)
        e = Integer( argv[4] );
    else
        std::cin >> e;

    IntRSADom<GivRandom> IR(n,e);

    std::ifstream TXT(argv[1]);
    if (!TXT) { std::cerr << "Error opening input file: " << argv[1] << std::endl; return -1; }
    std::ofstream OUT(argv[2]);
    if (!OUT) { std::cerr << "Error opening output file: " << argv[2] << std::endl; return -1; }
    tim.start();
    IR.encipher( OUT, TXT );
    OUT.close();
    TXT.close();
    tim.stop();

    std::cerr << tim << std::endl;


    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
