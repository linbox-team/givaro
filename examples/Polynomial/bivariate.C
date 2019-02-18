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
#include <givaro/modular-integer.h>
#include <givaro/givpoly1.h>
#include <givaro/givpower.h>
#include <givaro/givquotientdomain.h>
#include <string>
#include <sstream>

using namespace Givaro;


int main(int argc, char ** argv) {

    Modular<Integer> Zp( "1234567891234567919" );  // integers modulo 1234567891234567919
    typedef Modular<Integer> Field;
    typedef Poly1Dom< Field, Dense> Polynomials;
    typedef Poly1Dom< Polynomials, Dense> Bivariates;


    // Polynomials over Z13, with X as indeterminate
    Polynomials PZp( Zp, Indeter("X") );
    Bivariates  PPZp( PZp, Indeter("Y") );

    Field::Element tmp;
    Polynomials::Element P, R, S;
    PZp.init(P, Degree(1), 5 ); // 5X, PZp.init will call Zp.init on 5
    PZp.addin(P, Zp.init(tmp,7) ); // 5X+7

    std::istringstream stream(std::string("3 1 2 3 4"));
    PZp.read( stream, S );
    PZp.write(std::cout << "S: ", S ) << std::endl;

    PZp.assign(R, Degree(3), Zp.init(tmp,11) ); // 11X^3, PZp.assign will call Zp.assign on tmp
    PZp.addin(R, Zp.init(tmp,13) ); // 11X^3+13

    Bivariates::Element Q;
    PPZp.init(Q, Degree(2));
    PZp.assign(Q[0], R);
    PZp.assign(Q[2], P);// (5X+7)Y^2 + 11X^3+13

    PPZp.write(std::cout << "Q: ", Q) << std::endl;


    typedef QuotientDom<Bivariates> BivMods;
    BivMods QD(PPZp, Q);
    QD.write(std::cout << "Quotient: ") << std::endl;

    BivMods::Element Res, G; QD.init(Res); QD.init(G);

    PZp.assign(P, Degree(1), 1 ); 	// X
    PPZp.init(G, Degree(1), 1); 	// Y
    PPZp.addin(G, P);				// Y+X
    QD.write(std::cout << "Y+X: ", G) << std::endl;

    long l = 5;
    dom_power(Res, G, l, QD); 		// G^l mod Q

    QD.write(std::cout << "(Y+X)^" << l << ": ", Res) << std::endl;

    return 0;

}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
