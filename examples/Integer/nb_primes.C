// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

/*! @file examples/Integer/nb_primes.C
 * @ingroup examples
 * @ingroup integers
 * @example examples/Integer/nb_primes.C
 * @brief NO DOC
 */
#include <iostream>
#include "givaro/givintprime.h"
#include "givaro/givtimer.h"


using namespace Givaro;



inline Integer GIVMAX(const Integer& a, const Integer& b) { return (a<b?b:a); }

int main(int argc, char** argv)
{
    IntPrimeDom IPD;
    Integer m, tp = 2, ttp=1, np;
    if (argc > 1)
        np = Integer(argv[1]);
    else
        std::cin >> np;
    if (argc > 2)
        m = GIVMAX(1,Integer(argv[2]));
    else
        m = 2;
    Integer nf = m;
    unsigned long long nb = (IPD.isprime(m)?1:0);
    Timer tim; tim.clear(); tim.start();
    for (;m < np; tp *= 2) {
        std::cout << nb << " primes between " << nf << " and " << m << "\t\t == nextprime(" << (ttp+nf-1) << ")" << std::endl; ttp = tp;
        for (;(m < np) && (m < (nf+tp)); ++nb)
            IPD.nextprimein(m);
    }
    tim.stop();
    std::cout << (m>np?nb-1:nb) << " primes between " << nf << " and " << np << std::endl << tim << std::endl;
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
