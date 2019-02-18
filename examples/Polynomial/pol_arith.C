// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
/*! @file examples/Polynomial/pol_arith.C
 * @ingroup examples
 * @ingroup polynomials
 * @example examples/Polynomial/pol_arith.C
 * @brief NO DOC
 */

#include <iostream>
#include <givaro/gfq.h>
#include <givaro/givpoly1.h>

using namespace Givaro;


int main(int argc, char ** argv) {

    {
        GFqDom<int> Z13( 13, 1 );  // integers modulo 13
        GFqDom<int>::Element tmp;

        // Polynomials over Z13, with X as indeterminate
        Poly1Dom< GFqDom<int>, Dense > DP13( Z13, Indeter("X") );
        Poly1Dom< GFqDom<int>, Dense>::Element P, Q, R, monomial;

        // DP13.read( std::cin, P); // would read P as a succession of integers :
        // deg leadcoeff (lead-1)coeff ... unitcoeff

        DP13.init(P, {5,-33,12}); // P is now 5-33*X+12*X^2
        DP13.write( std::cout << "P: " , P )<< std::endl;

        DP13.assign( Q, Z13.init(tmp,6) ); // Q is degree 0 polynomial : 6 modulo 13
        DP13.write( std::cout << "Q: " , Q )<< std::endl;
        DP13.init( monomial, Degree(4), 3U);
        DP13.write( std::cout << "m: " , monomial )<< std::endl;
        DP13.addin( Q, monomial) ;
        DP13.write( std::cout << "Q: " , Q )<< std::endl;
        DP13.init( monomial, Degree(1), 75U);
        DP13.write( std::cout << "m: " , monomial )<< std::endl;
        DP13.addin( Q, monomial) ;
        DP13.write( std::cout << "Q: " , Q )<< std::endl;
        DP13.init( monomial, Degree(3), 45U);
        DP13.write( std::cout << "m: " , monomial )<< std::endl;
        DP13.subin( Q, monomial) ;
        DP13.write( std::cout << "Q: " , Q )<< std::endl;
        // Q is now 3*X^4+75*X-45*X^3+6

        DP13.mul ( R, P, Q); // R = P*Q;

        DP13.write( DP13.write(
                               std::cout << "(" , P ) << ") * (", Q) << ")";

        DP13.write(std::cout << " = " , R) << std::endl;

        DP13.gcd ( R, P, Q); //

        DP13.write( DP13.write( DP13.write(
                                           std::cout << "gcd(", P ) << ",", Q) << ") = ", R) << std::endl;

        DP13.lcm ( R, P, Q); //
        DP13.write( DP13.write( DP13.write(
                                           std::cout << "lcm(", P ) << ",", Q) << ") = ", R) << std::endl;
        DP13.lcm ( R, Q, P); //
        DP13.write( DP13.write( DP13.write(
                                           std::cout << "lcm(", Q ) << ",", P) << ") = ", R) << std::endl;

    }
    return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
