/*  ruint_arazi.cpp - Arazi & Qi for RecInt library test file

    Return value.
        0    No error
        != 0 Bad result for an operation

    The following constants have to be defined.
        STD_RECINT_SIZE     size of recint (> 6)
        LOOPS               number of loops
*/ 

#include <recint/ruint.h>
#include <recint/rmgmodule.h> /* arazi_qi() */

#if not defined(LOOPS)
#define LOOPS 1000
#endif

using namespace RecInt;

int main(void)
{
    ruint<STD_RECINT_SIZE> p, pinv, paz;
    ruint<STD_RECINT_SIZE+1> P, P1, R;
    
    // Random
    RecInt::srand(limb(time(NULL)));
    
    for (UDItype l = 0; l < LOOPS; l++) {
        rand(p); if (p % 2 == 0) p++;

        // P = p and R = 2^(2^K)
        copy(P.Low, p);
        R.High = 1;
        
        // P1 = inv(P) mod R
        inv_mod(P1, P, R);
        copy(pinv, P1.Low);
        
        // paz = inv(p) mod R
        arazi_qi(paz, p);

        if (paz != pinv) return 1;
    }
    
    return 0; 
}

