// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

/*! @file examples/Integer/primitiveroot.C
 * @ingroup examples
 * @ingroup integers
 * @example examples/Integer/primitiveroot.C
 * @brief NO DOC
 */
#include <iostream>
#define GIVARO_LENSTRA
#include <givaro/givintnumtheo.h>
#include <givaro/givtimer.h>


#ifndef TIMING
#define TIMING 1
#endif

using namespace Givaro;



int main(int argc, char** argv)
{
    IntNumTheoDom<> IP;
#ifndef __GIVARO_GMP_NO_CXX
    IP.seeding();
#endif
    IntNumTheoDom<>::Element a,pr;
    if (argc > 1) a = IntNumTheoDom<>::Element(argv[1]); else std::cin >> a;

    uint64_t runs;
    Timer tim; tim.clear();
    if (IP.isprime(a)) {
        Integer phin; IP.sub(phin,a,IP.one);
        std::vector<Integer> Lf;
        IP.write(std::cout << "Totient : ", Lf,phin) << std::endl;
        tim.start();
        for(uint64_t i = 0; i < TIMING; ++i)
            IP.prim_root_of_prime(pr, a);
        tim.stop();
        IP.write( std::cout << "Deterministic   : ", pr ) << std::endl;
        std::cerr << tim << std::endl;
    }
    tim.start();
    for(uint64_t i = 0; i < TIMING; ++i)
        IP.prim_root(pr, runs, a);
    tim.stop();
    IP.write( std::cout << "Random : ", pr ) << std::endl;
    std::cerr << tim << " (" << runs << " runs)" << std::endl;
    tim.start();
    for(uint64_t i = 0; i < TIMING; ++i)
        IP.lowest_prim_root(pr, a);
    tim.stop();
    IP.write( std::cout << "Lowest : ", pr ) << std::endl;
    std::cerr << tim << std::endl;

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
