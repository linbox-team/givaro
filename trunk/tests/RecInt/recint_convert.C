/* test_rand.cpp - Randomness of RecInt library test file

Return value.
    0    No error
    != 0 Bad result for an operation

The following constants have to be defined.
    STD_RECINT_SIZE     size of recint (> 5)
    LOOPS           number of loops of randomized tests
*/

#include <gmpxx.h>

#include <recint/recint.h>
#include <recint/ruconvert.h>

using namespace RecInt;

int main(void)
{
    ruint<STD_RECINT_SIZE> a, ca, p;
    rmint<STD_RECINT_SIZE> b, cb;
    mpz_class ga, gb;
      
    // Init.
    RecInt::srand(time(NULL));
    do { rand(p); } while(p % 2 == 0);
    b.init_module(p);
    
    // With ruint
    rand(a);
    ruint_to_mpz(ga, a);
    mpz_to_ruint(ca, ga);
    if (a != ca) return 1;
    
    fill_with_1(a);
    ruint_to_mpz(ga, a);
    mpz_to_ruint(ca, ga);
    if (a != ca) return 2;    

    // With rmint
    rand(b);
    rmint_to_mpz(gb, b);
    mpz_to_rmint(cb, gb);
    if (b != cb) return 8;
    
    b = p - 1;
    rmint_to_mpz(gb, b);
    mpz_to_rmint(cb, gb);
    if (b != cb) return 9;

    return 0; 
}

