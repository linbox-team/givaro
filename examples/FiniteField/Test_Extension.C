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
#include "givaro/extension.h"

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
template void FaireEssai< GFqDom<int64_t> >(const GFqDom<int64_t> & F) ;


int main (int argc, char * * argv) {

    uint64_t q = (argc>1?(uint64_t)atoi(argv[1]):13);
    uint64_t expo = (argc>2?(uint64_t)atoi(argv[2]):8);
    /*
       GFqDom<int64_t> Toto(q,1);
       Toto.write( std::cout << "This is ") << std::endl ;
       FaireEssai( Toto );
       */
    std::cerr << "Exponent max for zech logs with characteristic " << q << " : " << FF_EXPONENT_MAX(q,expo) << std::endl;
    std::cerr << "Sub-Exponent max for zech logs " << q << "^" << expo << " : " << FF_SUBEXPONENT_MAX(q,expo) << std::endl;
    std::cerr << "NEED polynomial representation : " << NEED_POLYNOMIAL_REPRESENTATION(q,expo) << std::endl;
    if ( NEED_POLYNOMIAL_REPRESENTATION(q,expo) )
        FaireEssai( Extension<>(q, expo) );
    else
        FaireEssai( GFqDom<int64_t>(q, expo) );

    FaireEssai( EXTENSION(q, expo) );
    return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
