#include <iostream>
#include <cstdlib>

#include <recint/recint.h>
#include <givaro/givtimer.h>

#define STD_RECINT_SIZE 9
#define LOOPS 500

#define ALEA_MAX  64
#define ALEA_MASK 63

using namespace RecInt;

int main(void)
{
    rmint<STD_RECINT_SIZE> m[ALEA_MAX];
    ruint<STD_RECINT_SIZE> u[ALEA_MAX], module;
    Givaro::Timer tim;
    
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
    
	// -----------
	// Standard output for benchmark - Alexis Breust 2014/12/11
	std::cout << "Time: " << tim.usertime()
			  << " Gflops: " << "Irrelevant" << std::endl;
    
    return 0;
}

