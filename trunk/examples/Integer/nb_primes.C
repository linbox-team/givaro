#include <iostream>
#include "givaro/givintprime.h"
#include "givaro/givtimer.h"


inline Integer GIVMAX(const Integer& a, const Integer& b) { return (a<b?b:a); }

int main(int argc, char** argv)
{
    IntPrimeDom IPD;
    Integer m, tp = 2, ttp=1, np;
    if (argc > 1)
        np = Integer(argv[1]);
    else
        std::cin >> np;
    if (argc > 2)
	m = GIVMAX(1,Integer(argv[2]));
    else
	m = 2;
    Integer nf = m;
    unsigned long long nb = (IPD.isprime(m)?1:0); 
    Timer tim; tim.clear(); tim.start();
    for (;m < np; tp *= 2) {
        std::cout << nb << " primes between " << nf << " and " << m << "\t\t == nextprime(" << (ttp+nf-1) << ")" << std::endl; ttp = tp;
        for (;(m < np) && (m < (nf+tp)); ++nb) 
	            IPD.nextprimein(m);
    }    
    tim.stop();
    std::cout << (m>np?nb-1:nb) << " primes between " << nf << " and " << np << std::endl << tim << std::endl;
    return 0;
}

