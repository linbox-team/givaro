#include <iostream>
#include <cstdlib>

#include <recint/recint.h>
#include <givaro/givtimer.h>
#include <givaro/givinteger.h>

#if not defined(STD_RECINT_SIZE)
#define STD_RECINT_SIZE 8
#endif

#if not defined(LOOPS)
#define LOOPS 1000000
#endif

#define ALEA_MAX  64
#define ALEA_MASK 63

using namespace RecInt;

int main(int argc, char ** argv)
{
    size_t nbloops = static_cast<size_t>((argc > 1)? atoi(argv[1]) : LOOPS);
    //    std::cerr << "nbloops: " << nbloops << std::endl;

    ruint<STD_RECINT_SIZE> a[ALEA_MAX], d[ALEA_MAX];
    mpz_class b[ALEA_MAX], c[ALEA_MAX];
    Givaro::Timer tim,gmp;


    // Randomness
    for (unsigned int i = 0; i < ALEA_MAX; i++) {
        rand(a[i]);
        ruint_to_mpz(b[i],a[i]);
    }

    // Main loop
    tim.clear(); tim.start();
    for (unsigned int l = 0; l < nbloops; l++) {
        d[l & ALEA_MASK] = a[l & ALEA_MASK];
        d[l & ALEA_MASK] *= a[l & ALEA_MASK];
    }
    tim.stop();

    // Main loop
    gmp.clear(); gmp.start();
    for (unsigned int l = 0; l < nbloops; l++) {
        c[l & ALEA_MASK] = b[l & ALEA_MASK];
        c[l & ALEA_MASK] *= b[l & ALEA_MASK];
    }
    gmp.stop();

    ruint<STD_RECINT_SIZE> module; rand(module);

    // -----------
    // Standard output for benchmark - Alexis Breust 2014/12/11
    std::cout
    << "SIZE: " << STD_RECINT_SIZE
    << " Time: " << tim.usertime() << ' ' << gmp.usertime()
    << " Mflops: " << std::scientific << (double(nbloops))/tim.usertime()/1000.0/1000.0 << ' ' << (double(nbloops))/gmp.usertime()/1000.0/1000.0
    << ' ' << a[(int)(module)& ALEA_MASK] << ' ' << b[(int)(module)& ALEA_MASK] << std::endl ;

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
