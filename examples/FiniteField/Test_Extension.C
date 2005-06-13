using namespace std;
#include "givaro/givpoly1.h"
#include "givaro/givextension.h"


template<class FField>
void FaireEssai(const FField & F) {

  F.write( cout << "Working in : " ) << endl;

  typename FField::element a, b, r;

  std::cout << "Enter an element of this field: ";  F.read( std::cin , a ); 
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
    
    cerr << "Exponent max for zech logs " << q << " : " << FF_EXPONENT_MAX(q) << endl;
    cerr << "NEED polynomial representation : " << NEED_POLYNOMIAL_REPRESENTATION(q,expo) << endl;

    FaireEssai( Toto );
    if ( NEED_POLYNOMIAL_REPRESENTATION(q,expo) )
        FaireEssai( Extension<>(q, expo) );
    else
        FaireEssai( GFqDom<long>(q, expo) );
       
      
    return 0;
}
