#include <iostream>
#include <cstdlib>
#include <gmpxx.h>

#include <recint/recint.h>
#include <givaro/givtimer.h>

#if not defined(STD_RECINT_SIZE)
#define STD_RECINT_SIZE 9
#endif

#if not defined(LOOPS)
#define LOOPS 1000000
#endif

#define ALEA_MAX  64
#define ALEA_MASK 63

using namespace RecInt;

int main(void)
{
    mpz_class ga[ALEA_MAX], gmodule;
    ruint<STD_RECINT_SIZE> module;
    rmint<STD_RECINT_SIZE> a;
    Givaro::Timer tim;
    
    std::cout << "Size of numbers: 2^(2^" << STD_RECINT_SIZE << ")" << std::endl;
    
    // For montgomery algorithm, the module must be odd
    RecInt::srand(42);
    rand(module);
    if (module % 2 == 0) module++;
    rmint<STD_RECINT_SIZE>::init_module(module);
    ruint_to_mpz(gmodule, module);
    
    // Randomness
    for (unsigned int i = 0; i < ALEA_MAX; i++) {
        rand(a);
        rmint_to_mpz(ga[i], a);
    }
    
    // Main loop
	tim.clear(); tim.start();
    for (unsigned int l = 0; l < LOOPS; l++) {
        mpz_mul(ga[l & ALEA_MASK].get_mpz_t(), ga[l & ALEA_MASK].get_mpz_t(), ga[(l+1) & ALEA_MASK].get_mpz_t());
        mpz_mod(ga[l & ALEA_MASK].get_mpz_t(), ga[l & ALEA_MASK].get_mpz_t(), gmodule.get_mpz_t());
    }
    tim.stop();
    
    std::cout << "RecInt: " << tim.usertime() << std::endl;
    std::cout << "Result: " << std::hex << ga[0] << std::endl;
    
    return 0;
}

