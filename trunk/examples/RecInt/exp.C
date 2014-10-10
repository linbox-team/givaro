#include <iostream>
#include <cstdlib>

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
    rmint<STD_RECINT_SIZE> m[ALEA_MAX];
    ruint<STD_RECINT_SIZE> u[ALEA_MAX], module;
    Givaro::Timer tim;
    
    std::cout << "Size of numbers: 2^(2^" << STD_RECINT_SIZE << ")" << std::endl;
    
    // For montgomery algorithm, the module must be odd
    RecInt::srand(42);
    rand(module);
    if (module % 2 == 0) module++;
    rmint<STD_RECINT_SIZE>::init_module(module);
    
    // Randomness
    for (unsigned int i = 0; i < ALEA_MAX; i++) {
        rand(m[i]);
        rand(u[i]);
    }
    
    // Main loop
	tim.clear(); tim.start();
    for (unsigned int l = 0; l < LOOPS; l++) {
        exp(m[l & ALEA_MASK],     m[(l+2) & ALEA_MASK], u[l & ALEA_MASK]);
        exp(m[(l+2) & ALEA_MASK], m[(l+1) & ALEA_MASK], u[l & ALEA_MASK]);
    }
    tim.stop();
    
    std::cout << "RecInt: " << tim.usertime() << std::endl;
    std::cout << "Result: " << std::hex << m[0] << std::endl;
    
    return 0;
}

