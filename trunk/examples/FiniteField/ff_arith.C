#include <iostream>
#include <givaro/givzpz.h>
#include <givaro/givgfq.h>

int main(int argc, char ** argv) {

 {
    ZpzDom<Std32> Z13(13);   // modulo 13 over 32 bits
    ZpzDom<Std32>::element a, b, c;
    Z13.init(a, 7);
    Z13.init(b, -29.0);
    
    Z13.add(c, a, b); // c = a+b

    std::cerr << "Within "; 
    Z13.write( std::cerr );
    std::cerr << " : " << std::flush;

     // Separate output writing
    Z13.write( std::cout, a) << " + " << std::flush;  
    Z13.write( std::cout, b) << " = " << std::flush; 
    Z13.write( std::cerr, c) << std::endl;
 }


 {
    int Mod = 3; int exponent = 4;
    GFqDom<int> GF81( Mod, exponent );  // finite field with 3^4 elements
    GFqDom<int>::element a, b, c;

    GF81.init(a, 7);    // 7 modulo   3 
    GF81.init(b, -28);  // -28 modulo 3

    GF81.mul(c, a, b);  // c = a*b

    GFqDom<int>::element d; 

    GF81.axpy(d, c, a, b); // d = c + a*b;

    GF81.write( std::cerr << "Within " ) << " : " << std::flush;

    // Writing all outputs in a single command line
    GF81.write( GF81.write( GF81.write( GF81.write(
       std::cout, c) << " + ", a) << " * ", b) << " = ", d) << std::endl;

 }   

 {
    ZpzDom<Integer> Z13(13);   // modulo 13 over arbitrary size
    ZpzDom<Integer>::element a, b, c;
    Z13.init(a, 7);
    Z13.init(b, -29);
    

    std::cerr << "Within "; 
    Z13.write( std::cerr );
    std::cerr << " : " << std::flush;

    Z13.write( Z13.write( Z13.write( std::cout, a) << " + ",  b) << " = ", 
    	Z13.add(c, a, b)
	     ) << std::endl;
    Z13.write( Z13.write( Z13.write( std::cout, a) << " - ",  b) << " = ", 
    	Z13.sub(c, a, b)
	     ) << std::endl;
    Z13.write( Z13.write( Z13.write( std::cout, a) << " * ",  b) << " = ", 
    	Z13.mul(c, a, b)
	     ) << std::endl;
    Z13.write( Z13.write( Z13.write( std::cout, a) << " / ",  b) << " = ", 
    	Z13.div(c, a, b)
	     ) << std::endl;
 }

    return 0;
}

