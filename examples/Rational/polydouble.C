// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

/*! @file examples/Rational/polydouble.C
 * @ingroup examples
 * @ingroup rationals
 * @example examples/Rational/polydouble.C
 * @brief NO DOC
 */
#include <iostream>
#include <stdlib.h>
#include <givaro/qfield.h>
#include <givaro/givpoly1.h>
#include <givaro/givtimer.h>
#include <givaro/givinit.h>         // Givaro initialization
#include <givaro/givprint.h>        // Givaro print utils

using namespace Givaro;

typedef Poly1Dom< QField<Rational>, Dense>::Element RatPoly;
typedef std::vector<double> DoublePoly;

int main(int argc, char** argv)
{
    srandom((unsigned int)BaseTimer::seed() );


    Integer f,m,k;

    QField<Rational> Q;
    Poly1Dom< QField<Rational>, Dense> PolQ(Q);


    DoublePoly D;
    RatPoly R;


    size_t n = (argc>1?(size_t)atoi(argv[1]):10);

    for(size_t i=0;i<n;++i)
        D.push_back( (double(random()) / RAND_MAX) );


    R.resize( D.size() );

    RatPoly::iterator it=R.begin();
    DoublePoly::const_iterator dit = D.begin();
    for( ; dit != D.end(); ++dit, ++it)
        *it =*dit;

    std::cout << "Double Poly : " << D << std::endl;
    std::cout << "REDUCED Rational " << R << std::endl;
    std::cout << "Approximations : ";

    RatPoly::const_iterator cit=R.begin();
    dit = D.begin();
    for( ; dit != D.end(); ++dit, ++cit)
        std::cout << std::endl << *cit << " is " << ((double)*cit) << " by " << ( (double)*cit -*dit ) << ' ';
    std::cout << std::endl;



    Rational::SetNoReduce();

    it=R.begin();
    dit = D.begin();
    for( ; dit != D.end(); ++dit, ++it)
        *it =*dit;

    std::cout << "Double Poly : " << D << std::endl;
    std::cout << "Unreduced Rational " << R << std::endl;
    std::cout << "Approximations : ";

    cit=R.begin();
    dit = D.begin();
    for( ; dit != D.end(); ++dit, ++cit)
        std::cout << std::endl << *cit << " is " << ((double)*cit) << " by " << ( (double)*cit -*dit ) << ' ';
    std::cout << std::endl;



    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
