#include <iostream>
#include <givaro/givzpz64std.h>
#include <givaro/givzpz.h>
#include <givaro/givgfq.h>
#include <givaro/givmontg32.h>


template<class Field>
void TestField(const Field& F) {
    std::cerr << "Within " ; 
    F.write(std::cerr );
    std::cerr  << " : " << std::flush;

    typename Field::Element a, b, c, d;
    F.init(a, 7);
    F.init(b, -29.3);

    F.init(c);            // empty constructor
    F.init(d);            // empty constructor
    
    F.add(c, a, b);       // c = a+b
    
     // Separate output writing
    F.write( std::cout, a) << " + " << std::flush;  
    F.write( std::cout, b) << " = " << std::flush; 
    F.write( std::cerr, c) << std::endl;
    

    F.mul(c, a, b);     // c = a*b
    F.axpy(d, a, b, c); // d = a*b + c;

    // Writing all outputs in a single command line
    F.write( std::cerr << "Within " ) << " : " << std::flush;
    F.write( F.write( F.write( F.write(
        std::cout, c) << " + ", a) << " * ", b) << " = ", d) << std::endl;


        // Four operations
    F.write( F.write( F.write( std::cout, a) << " += ",  b) << " is ", 
             F.addin(a, b)
	     ) << "   ;   ";
    F.write( F.write( F.write( std::cout, a) << " -= ",  b) << " is ", 
             F.subin(a, b)
	     ) << "   ;   ";
    F.write( F.write( F.write( std::cout, a) << " *= ",  b) << " is ", 
             F.mulin(a, b)
	     ) << "   ;   ";
    F.write( F.write( F.write( std::cout, a) << " /= ",  b) << " is ", 
             F.divin(a, b)
	     ) << std::endl;
}



int main(int argc, char ** argv) {

        // modulo 13 over 16 bits
    ZpzDom<Std16> C13(13); TestField( C13 );
    
        // modulo 13 over 32 bits
    ZpzDom<Std32> Z13(13); TestField( Z13 );

        // modulo 13 over unsigned 32 bits
    ZpzDom<Unsigned32> U13(13); TestField( U13 );

#ifdef __USE_Givaro_SIXTYFOUR__
        // modulo 13 over 64 bits
    ZpzDom<Std64> LL13(13UL); TestField( LL13 );
#endif

        // modulo 13 fully tabulated
    ZpzDom<Log16> L13(13); TestField( L13 );

        // modulo 13 over 32 bits with Montgomery reduction
    Montgomery<Std32> M13(13); TestField( M13 );
    
        // modulo 13 with primitive root representation
    GFqDom<int> GF13( 13 ); TestField( GF13 );
    
        // finite field with 5^4 elements
    GFqDom<int> GF81( 5, 4 ); TestField( GF81 );
    
        // modulo 13 over arbitrary size
    ZpzDom<Integer> IntZ13(13); TestField( IntZ13 );

    

    return 0;
}

