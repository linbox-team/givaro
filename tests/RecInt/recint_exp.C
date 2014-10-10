/* rmint_inv.cpp - Inverse in modular arithmetic of RecInt library test file

Return value.
    0    No error
    != 0 Bad result for an operation

The following constants have to be defined.
    STD_RECINT_SIZE     size of recint (> 5)
    LOOPS           number of loops of randized tests
*/

#include <gmpxx.h>

#include <recint/recint.h>

#if not defined(LOOPS)
#define LOOPS 1000
#endif

using namespace RecInt;

int main(void)
{
    rmint<STD_RECINT_SIZE> a, b;
    ruint<STD_RECINT_SIZE> c, au, bu;
    mpz_class ga, gb, gc, gp, gcmp;
    USItype r;
      
    // Init. size = p
    RecInt::srand(time(NULL));
    ruint<STD_RECINT_SIZE> p;

    // Loop
    for (UDItype l = 1; l < LOOPS; l++) {
        do { rand(p); } while ((p % 2) == 0);
        a.init_module(p);
        ruint_to_mpz(gp, p);
        
        //------- Exp ---------
        
        // With USItype
        rand(b); r = rand();
        rmint_to_mpz(gb, b);
        exp(a, b, r);
        mpz_powm_ui(ga.get_mpz_t(), gb.get_mpz_t(), r, gp.get_mpz_t());
        rmint_to_mpz(gcmp, a);
        if (gcmp != ga) return 1;
        
        // With ruint
        rand(c); r = rand();
        ruint_to_mpz(gc, c);
        exp(a, b, c);
        mpz_powm(ga.get_mpz_t(), gb.get_mpz_t(), gc.get_mpz_t(), gp.get_mpz_t());
        rmint_to_mpz(gcmp, a);
        if (gcmp != ga) return 2;
        
        // Exp mod (all in ruint)
        rand(au); rand(bu);
        ruint_to_mpz(ga, au);
        ruint_to_mpz(gb, bu);
        exp_mod(c, au, bu, p);
        mpz_powm(gc.get_mpz_t(), ga.get_mpz_t(), gb.get_mpz_t(), gp.get_mpz_t());
        ruint_to_mpz(gcmp, c);
        if (gcmp != gc) return 3;
    }
        
    return 0; 
}

