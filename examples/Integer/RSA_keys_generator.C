// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

/*! @file examples/Integer/RSA_keys_generator.C
 * @ingroup examples
 * @ingroup integers
 * @example examples/Integer/RSA_keys_generator.C
 * @brief NO DOC
 */
#include <iostream>
#include "givaro/givintrsa.h"
#include "givaro/givrandom.h"
#include "givaro/givtimer.h"

using namespace Givaro;




int main(int argc, char** argv)
{
    Timer tim;
    tim.clear();

    // m will be p times q
    // Sizes are in number of int to compose the big integer
    long keysize;
    if (argc > 1)
        keysize = atoi( argv[1] );
    else
        std::cin >> keysize;


    GivRandom generator;
    Integer::seeding(generator.seed());


    tim.start();
    IntRSADom<GivRandom> IR(keysize,true,generator);
    tim.stop();

    std::cout << IR.getn() << " " << IR.gete() << " " << IR.getd()  << std::endl;
    std::cerr << tim << std::endl;

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
