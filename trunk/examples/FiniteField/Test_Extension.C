// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software. 
// see the COPYRIGHT file for more details.

using namespace std;
#include "givaro/givpoly1.h"
#include "givaro/givextension.h"


template<class FField>
void FaireEssai(const FField & F) {

  F.write( cout << "Working in : " ) << endl;

  typename FField::Element a, b, r;

  std::cout << "Enter an Element of this field: ";  F.read( std::cin , a ); 
//   F.init( a, "1+3*X+5*X^2" );
  F.init( b, (Integer)23 );

  F.add(r, a, b);

  F.write( F.write( F.write( F.write(std::cout, a) << " + " , b ) << " = " , r) << " with ") << endl ; 

}



int main (int argc, char * * argv) {
    
    unsigned long q = (argc>1?atoi(argv[1]):13);
    unsigned long expo = (argc>2?atoi(argv[2]):8);

    GFqDom<long> Toto(q,1);
    Toto.write( cout << "This is ") << endl ;
    FaireEssai( Toto );

    cerr << "Exponent max for zech logs " << q << "^" << expo << " : " << FF_EXPONENT_MAX(q,expo) << endl;
    cerr << "NEED polynomial representation : " << NEED_POLYNOMIAL_REPRESENTATION(q,expo) << endl;
    if ( NEED_POLYNOMIAL_REPRESENTATION(q,expo) )
        FaireEssai( Extension<>(q, expo) );
    else
        FaireEssai( GFqDom<long>(q, expo) );
       
      
    return 0;
}
