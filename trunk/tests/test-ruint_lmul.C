/* test_lmul.cpp - Arithmetic multiplications of RecInt generic test file

Return value.
    0    No error
    != 0 Bad result for an operation

The following constants have to be defined.
    STD_RECINT_SIZE     size of recint (> 5)
    LOOPS           number of loops of randized tests
*/

#include <gmpxx.h>
#include <recint/ruint.h>

#if not defined(LOOPS)
#define LOOPS 10000
#endif

using namespace RecInt;

int main(void)
{
    ruint<STD_RECINT_SIZE+1> z;
    ruint<STD_RECINT_SIZE> x, y, zh, zl;
    mpz_class gx, gy, gz, gcmp;
      
    // Init.
    RecInt::srand(limb(time(NULL)));
    
    // Loop
    for (UDItype l = 1; l < LOOPS; l++) {
        // RecInt rand
        rand(x);
        rand(y);
        rand(z);
        ruint_to_mpz(gx, x);
        ruint_to_mpz(gy, y);
        gz = gx * gy;
        
        // Naive method
        lmul_naive(z, x, y);
        ruint_to_mpz(gcmp, z);
        if (gcmp != gz) return 1;
        
        lmul_naive(zh, zl, x, y);
        z.High = zh; z.Low = zl;
        ruint_to_mpz(gcmp, z);
        if (gcmp != gz) return 1;
        
        // Karatsuba method
        lmul_kara(z, x, y);
        ruint_to_mpz(gcmp, z);
        if (gcmp != gz) return 2;
        
        lmul_kara(zh, zl, x, y);
        z.High = zh; z.Low = zl;
        ruint_to_mpz(gcmp, z);
        if (gcmp != gz) return 2;
        
    }

    return 0; 
}

