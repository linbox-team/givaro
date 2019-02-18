// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

/*! @file examples/Integer/probable_primroot.C
 * @ingroup examples
 * @ingroup integers
 * @example examples/Integer/probable_primroot.C
 * @brief NO DOC
 */
#include <iostream>
#include <givaro/givintnumtheo.h>
#include <givaro/givtimer.h>
#include <cmath>

using namespace Givaro;



//  Polynomial-time generation of primitive roots
//  L is number of loops of Pollard partial factorization of n-1
//  10,000,000 gives at least 1-2^{-40} probability of success
//  [Dubrois & Dumas, Industrial-strength primitive roots]
//  Returns the probable primitive root and the probability of error.

int main(int argc, char** argv)
{
    IntNumTheoDom<> IP;
#ifdef __GMP_PLUSPLUS__
    IP.seeding( (unsigned long)BaseTimer::seed() );
#endif

    double error;
    IntNumTheoDom<>::Element a,pr,g;
    if (argc > 1) a = IntNumTheoDom<>::Element(argv[1]); else std::cin >> a;
    bool comp ; if ( (comp=(! IP.isprime(a))) )  std::cerr << a << " is not prime, primitive root will have no sense and may loop forever ..." << std::endl;
    double epsilon = argc > 2 ? atof(argv[2]) : 0.0000001;

    Timer tim; tim.clear();
    tim.start();
    //======================================================================
    // Default is partial factorization up to factors of at least 12 digits.
    // with probability of error much less than 2^{-40}
    // IP.probable_prim_root(pr, error, a );


    //======================================================================
    // Choosing L to be O(log^2(p))
    // gives the best probability with O(log^4(p)) complexity
    // Probability of error is approximately O(1/log^4(p))
    // IP.probable_prim_root(pr, error, a, (unsigned long)power(logtwo(a),2));

    //======================================================================
    // Choosing L to be O( \sqrt(epsilon) )
    // gives probability of error at most epsilon
    // Newton-Raphson iteration is used for
    // 1-epsilon = (1+2/(p-1))*(1-1/B)^(ln( (p-1)/2 )/ln(B))
    // So that no factor less than B can be avoided
    // With Pollard's rho factorization, L is chosen to be sqrt(B)
    // Might not be polynomial if epsilon is too big
#define GIVARO_POLLARD
    IP.probable_prim_root(pr, error, a, epsilon );
    tim.stop();

    if (comp) std::cerr << IP.gcd(g,a,pr) << " is a factor of " << a << std::endl;
    IntegerDom().write( std::cout << "Prim root   : ", pr );
    if (error > 0) {
        std::cout << ", correct with probability at least : 1-" << error << std::endl;
        std::cerr << tim << std::endl;

        std::cerr << "Now checking primitivity, this may take some time (complete factorization of n-1) ...";

#define GIVARO_LENSTRA

        Timer verif; verif.clear(); verif.start();
        if ( IP.isorder(a-1, pr, a) ) {
            verif.stop();
            std::cerr << "... Pimitivity checked" << std::endl;
            std::cerr << verif << std::endl;
        }
        else {
            verif.stop();
            std::cerr << "... WARNING : FAILURE" << std::endl;
            std::cerr << verif << std::endl;
        }
    } else {
        std::cout << ", deterministically correct" << std::endl;
        std::cerr << tim << std::endl;
    }

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
