#include <iostream>
#include <cstdlib>

#include <recint/recint.h>
#include <givaro/givtimer.h>

#if not defined(STD_RECINT_SIZE)
#define STD_RECINT_SIZE 9
#endif

#if not defined(LOOPS)
#define LOOPS 10000
#endif

#define ALEA_MAX  64
#define ALEA_MASK 63

using namespace RecInt;

int main(void)
{
    rmint<STD_RECINT_SIZE> a[ALEA_MAX];
    ruint<STD_RECINT_SIZE> module;
    Givaro::Timer tim;
    
    // For montgomery algorithm, the module must be odd
    RecInt::srand(42);
    rand(module);
    if (module % 2 == 0) module++;
    rmint<STD_RECINT_SIZE>::init_module(module);
    
    // Randomness
    for (unsigned int i = 0; i < ALEA_MAX; i++)
        rand(a[i]);
    
    // Main loop
	tim.clear(); tim.start();
    for (unsigned int l = 0; l < LOOPS; l++) {
        mul(a[l & ALEA_MASK], a[l & ALEA_MASK], a[(l+1) & ALEA_MASK]);
    }
    tim.stop();
    
	// -----------
	// Standard output for benchmark - Alexis Breust 2014/12/11
	std::cout << "Time: " << tim.usertime()
			  << " Gflops: " << "Irrelevant" << std::endl;
    
    return 0;
}

