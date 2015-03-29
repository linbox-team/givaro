// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
/*! @file examples/FiniteField/Test_Extension.C
 * @ingroup examples
 * @ingroup finitefields
 * @example examples/FiniteField/Test_Extension.C
 * @brief NO DOC
 */

#include "givaro/givpoly1.h"
#include "givaro/givextension.h"

using namespace Givaro;




template<class FField> void FaireEssai(const FField & F) {

  F.write( std::cout << "Working in : " ) << std::endl;

  typename FField::Element a, b, r;

  std::cout << "Enter an Element of this field: ";  F.read( std::cin , a );
//   F.init( a, "1+3*X+5*X^2" );
  F.init( b, (Integer)23 );

  F.add(r, a, b);

  F.write( F.write( F.write( F.write(std::cout, a) << " + " , b ) << " = " , r) << " with ") << std::endl ;

}

template void FaireEssai< Extension<> >(const Extension<> & F) ;
template void FaireEssai< GFqDom<long> >(const GFqDom<long> & F) ;


int main (int argc, char * * argv) {

    unsigned long q = (argc>1?(unsigned long)atoi(argv[1]):13);
    unsigned long expo = (argc>2?(unsigned long)atoi(argv[2]):8);
/*
    GFqDom<long> Toto(q,1);
    Toto.write( std::cout << "This is ") << std::endl ;
    FaireEssai( Toto );
*/
    std::cerr << "Exponent max for zech logs with characteristic " << q << " : " << FF_EXPONENT_MAX(q,expo) << std::endl;
    std::cerr << "Sub-Exponent max for zech logs " << q << "^" << expo << " : " << FF_SUBEXPONENT_MAX(q,expo) << std::endl;
    std::cerr << "NEED polynomial representation : " << NEED_POLYNOMIAL_REPRESENTATION(q,expo) << std::endl;
    if ( NEED_POLYNOMIAL_REPRESENTATION(q,expo) )
        FaireEssai( Extension<>(q, expo) );
    else
        FaireEssai( GFqDom<long>(q, expo) );

    FaireEssai( EXTENSION(q, expo) );
    return 0;
}
