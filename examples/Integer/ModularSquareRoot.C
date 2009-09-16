#include <iostream>
#include <stdlib.h>
#include <givaro/givintsqrootmod.h>
#include <givaro/givtimer.h>

// Algorithm 3.34 (Square Root Mod p) of 
// Handbook of Applied Cryptography
// by Menezes, van Oorschot, Vanstone

int main(int argc, char** argv)
{
    Integer a(argv[1]), n(argv[2]);
    std::cerr << "n: " << n << std::endl;
    std::cerr << "a: " << a << std::endl;
    
    
    if ( legendre( a , n ) == -1 ){ 
        std::cerr << "Not a quadratic residue" << std::endl;
   }
   else {
        IntSqrtModDom<> ISM;

	Integer r;
	Timer chrono; chrono.start();
	ISM.sqrootmod(r,a,n);
	chrono.stop();	
	std::cout << r << std::endl;
        std::cerr << chrono << std::endl;
        
        std::cerr << "Check, " << r << "^2 mod " << n << " = " << ( (r*r)%n) << std::endl;	

        }
        

    
    return 0;
}
