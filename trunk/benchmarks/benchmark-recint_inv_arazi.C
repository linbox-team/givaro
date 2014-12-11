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

using namespace RecInt;

int main(void)
{
    ruint<STD_RECINT_SIZE> p, pinv;
    Givaro::Timer tim;
    
    // Random
    RecInt::srand(time(NULL));
    
	tim.clear(); tim.start();
    for (UDItype l = 0; l < LOOPS; l++) {
        rand(p); if (p % 2 == 0) p++;

        arazi_qi(pinv, p);
    }
    tim.stop();
    
	// -----------
	// Standard output for benchmark - Alexis Breust 2014/12/11
	std::cout << "Time: " << tim.usertime()
			  << " Gflops: " << "Irrelevant" << std::endl;
    
    return 0; 
}

