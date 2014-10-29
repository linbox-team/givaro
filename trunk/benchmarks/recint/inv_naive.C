/*  ruint_arazi.cpp - Arazi & Qi for RecInt library test file

    Return value.
        0    No error
        != 0 Bad result for an operation

    The following constants have to be defined.
        STD_RECINT_SIZE     size of recint (> 6)
        LOOPS               number of loops
*/

#include <recint/ruint.h>
#include <givaro/givtimer.h>

#if not defined(STD_RECINT_SIZE)
#define STD_RECINT_SIZE 9
#endif

#if not defined(LOOPS)
#define LOOPS 100000
#endif

using namespace RecInt;

int main(void)
{
    ruint<STD_RECINT_SIZE> p, pinv;
    ruint<STD_RECINT_SIZE+1> P, P1, R;
    Givaro::Timer tim;
    
    // Random
    RecInt::srand(time(NULL));
    
    R.High = 1;
	tim.clear(); tim.start();
    for (UDItype l = 0; l < LOOPS; l++) {
        rand(P.Low); if (P.Low % 2 == 0) P.Low++;

        // P1 = inv(P) mod R
        inv_mod(P1, P, R);
        copy(pinv, P1.Low);
    }
    tim.stop();
    
    std::cout << "Old: " << tim.usertime() << std::endl;
    std::cout << "Last inv: " << std::hex << pinv << std::endl;
    
    return 0; 
}

