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
    Polynomials::Element P, R; 
    PZp.init(P, Degree(1), Zp.init(tmp,5) ); // 5X
    PZp.addin(P, Zp.init(tmp,7) ); // 5X+7

    PZp.init(R, Degree(3), Zp.init(tmp,11) ); // 11X^3
    PZp.addin(R, Zp.init(tmp,13) ); // 11X^3+13

    Bivariates::Element Q;
    PPZp.init(Q, Degree(2)); 
    PZp.assign(Q[0], R);
    PZp.assign(Q[2], P);// (5X+7)Y^2 + 11X^3+13
    

    PPZp.write(std::cout << "Q: ", Q) << std::endl;
    

    return 0;
    
}
