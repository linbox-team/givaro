#include <iostream>
#include <cstdlib>
#include <gmpxx.h>

#include <recint/recint.h>
#include <givaro/givtimer.h>

#if not defined(STD_RECINT_SIZE)
#define STD_RECINT_SIZE 9
#endif

#if not defined(LOOPS)
#define LOOPS 500
#endif

#define ALEA_MAX  64
#define ALEA_MASK 63

using namespace RecInt;

int main(void)
{
    mpz_class gm[ALEA_MAX], gu[ALEA_MAX], gmodule;
    ruint<STD_RECINT_SIZE> u, module;
    rmint<STD_RECINT_SIZE> m;
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
        rand(m);
        rand(u);
        rmint_to_mpz(gm[i], m);
        ruint_to_mpz(gu[i], u);
    }
    
    // Main loop
	tim.clear(); tim.start();
    for (unsigned int l = 0; l < LOOPS; l++) {
        mpz_powm(gm[l & ALEA_MASK].get_mpz_t(), gm[(l+2) & ALEA_MASK].get_mpz_t(), gu[l & ALEA_MASK].get_mpz_t(), gmodule.get_mpz_t());
        mpz_powm(gm[(l+2) & ALEA_MASK].get_mpz_t(), gm[(l+1) & ALEA_MASK].get_mpz_t(), gu[l & ALEA_MASK].get_mpz_t(), gmodule.get_mpz_t());
    }
    tim.stop();
    
    std::cout << "GMP: " << tim.usertime() << std::endl;
    std::cout << "Result: " << std::hex << gm[0] << std::endl;
    
    return 0;
}

