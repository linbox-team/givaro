#include <iostream>
#include <stdlib.h>
#include <givaro/givintprime.h>
#include <givaro/givrandom.h>
#include <givaro/givtimer.h>

// Algorithm 3.34 (Square Root Mod p) of 
// Handbook of Applied Cryptography
// by Menezes, van Oorschot, Vanstone

int main(int argc, char** argv)
{
    Integer a(argv[1]), p(argv[2]);
    std::cerr << "p: " << p << std::endl;
    std::cerr << "a: " << a << std::endl;
    if ( legendre( a , p ) == -1 ) 
        std::cerr << "Not a quadratic residue" << std::endl;
    else {
        Timer chrono; chrono.start();

        Integer t(p-1);
        size_t s=0;
        for( ; ((t>>1)<<1) == t ; ++s) t>>=1;
    
        size_t l=(size_t)ceil(logtwo(p));

        Integer::seeding( BaseTimer::seed() );
        Integer b; 
        while ( legendre( Integer::random(b,l), p ) != -1 ) ;

        Integer inva; inv(inva, a, p);

        Integer c; powmod(c, b, t, p);
        Integer r; powmod(r, a, (t+1)>>1, p);
        
        for(size_t i=1; i<s; ++i) {
            Integer d(r);
            d *= d; d %= p;
            d *= inva; d %= p;
            powmod(d, d, 1<<(s-i-1), p);
            if (d == p-1) {
                r *= c; r %= p;
            }
            c *= c ; c %= p;
        }
        chrono.stop();

        std::cout << r << std::endl;
        std::cerr << chrono << std::endl;
        
        std::cerr << "Check, " << r << "^2 mod " << p << " = " << ( (r*r)%p) << std::endl;
    }
    
    return 0;
}
