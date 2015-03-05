/*  ruint_arazi.cpp - Arazi & Qi for RecInt library test file

    Return value.
        0    No error
        != 0 Bad result for an operation

    The following constants have to be defined.
        STD_RECINT_SIZE     size of recint (> 6)
        LOOPS               number of loops
*/

#include <recint/ruint.h>
#include <recint/rmintmg.h> /* arazi_qi() */
#include <givaro/givtimer.h>

#if not defined(STD_RECINT_SIZE)
#define STD_RECINT_SIZE 9
#endif

#if not defined(LOOPS)
#define LOOPS 10000
#endif

#define ALEA_MAX  64
#define ALEA_MASK 63

using namespace RecInt;

int main(int argc, char ** argv)
{
    size_t nbloops = (argc>1?atoi(argv[1]):LOOPS);
    
    Givaro::Timer tim;
    
    ruint<STD_RECINT_SIZE> a[ALEA_MAX];
    ruint<STD_RECINT_SIZE> pinv[ALEA_MAX];
    // Randomness
    for (unsigned int i = 0; i < ALEA_MAX; i++) {
        rand(a[i]);
        if (a[i] % 2 == 0) ++a[i];
    }
    
 

    // Random
    RecInt::srand(time(NULL));
    
	tim.clear(); tim.start();
    for (UDItype l = 0; l < nbloops; l++) {
        arazi_qi(pinv[l&ALEA_MASK], a[l&ALEA_MASK]);
    }
    tim.stop();
    
	// -----------
	// Standard output for benchmark - Alexis Breust 2014/12/11
	std::cout << "Time: " << tim.usertime()
			  << " Gflops: " << std::scientific << (double(nbloops))/tim.usertime()/1000.0/1000.0/1000.0 << ' ' << a[(int)(rand(a[0]))& ALEA_MASK] << std::endl ;
    
    return 0; 
}

