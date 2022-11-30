// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

/*! @file examples/Integer/ifactor.C
 * @ingroup examples
 * @ingroup integers
 * @example examples/Integer/ifactor.C
 * @brief NO DOC
 */
#include <iostream>

#include <givaro/givinit.h>
#include <givaro/givintfactor.h>
#include <givaro/givtimer.h>

using namespace Givaro;

int main(int argc, char** argv)
{
    IntFactorDom<> IP;
    Integer m;
    if (argc > 1)
        m = Integer(argv[1]);
    else
        std::cin >> m;

    const size_t nbloops(argc>2?atoi(argv[2]):50000);

    if (IP.islt(m,0) ) {
        std::cout << "-";
        IP.negin(m);
    }

    if (IP.islt(m,4))
        IP.write(std::cout,m) << std::endl;
    else {
        std::vector<Integer> Moduli;
        std::vector<size_t> exponents;

        Timer tim; tim.clear(); tim.start();
        bool success = IP.set(Moduli, exponents, m, nbloops);
        tim.stop();

        for(size_t i=0; i<Moduli.size(); ++i)
            std::cout << Moduli[i] << '^' << exponents[i] << '*';
        std::cout << '1' << std::endl;

        std::clog << (success?"Complete ":"Incomplete ") << tim << std::endl;
    }
    return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
