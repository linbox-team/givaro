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
#include <givaro/givpoly1.h>
#include <givaro/givrational.h>
#include <givaro/givtimer.h>
#include <givaro/givinit.h>         // Givaro initialization

using namespace Givaro;

typedef Poly1Dom< RationalDom, Dense>::Element RatPoly;
typedef std::vector<double> DoublePoly;

std::ostream& operator<< (std::ostream& o, const RatPoly& v) {
    o << "Poly (s [v_1, ..,  v_s]) : ";
    o << v.size() << " [";
    for(size_t i=0; i<v.size(); ++i) {
        o << ' ' << v[i];
    }
    return o << ']';
}

std::ostream& operator<< (std::ostream& o, const DoublePoly& v) {
    o << '[';
    for(size_t i=0; i<v.size(); ++i) {
        o << ' ' << v[i];
    }
    return o << ']';
}

int main(int argc, char** argv)
{
    srandom(BaseTimer::seed() );


    Integer f,m,k;

    RationalDom Q;
    Poly1Dom< RationalDom, Dense> PolQ(Q);


    DoublePoly D;
    RatPoly R;


    size_t n = (argc>1?atoi(argv[1]):10);

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

