#include <iostream>
#include <cstdlib>

#include <recint/recint.h>
#include <givaro/givtimer.h>

#if not defined(STD_RECINT_SIZE)
#define STD_RECINT_SIZE 8
#endif

#if not defined(LOOPS)
#define LOOPS 1000
#endif

#define ALEA_MAX  64
#define ALEA_MASK 63

using namespace RecInt;

int main(int argc, char ** argv)
{
    size_t nbloops = static_cast<size_t>((argc > 1)? atoi(argv[1]) : LOOPS);

    rmint<STD_RECINT_SIZE> m[ALEA_MAX];
    ruint<STD_RECINT_SIZE> u[ALEA_MAX], module;
    mpz_class b[ALEA_MAX], c[ALEA_MAX], gmod;
    Givaro::Timer tim, gmp;

    // For montgomery algorithm, the module must be odd
    RecInt::srand(42);
    rand(module);
    if (module % 2 == 0) module++;
    rmint<STD_RECINT_SIZE>::init_module(module);
    ruint_to_mpz(gmod, module);

    // Randomness
    for (unsigned int i = 0; i < ALEA_MAX; i++) {
        rand(m[i]); ruint_to_mpz(b[i],m[i].Value);
        rand(u[i]); ruint_to_mpz(c[i],u[i]);
    }

    // Main loop
    tim.clear(); tim.start();
    for (unsigned int l = 0; l < nbloops; l++) {
        exp(m[l & ALEA_MASK],     m[(l+2) & ALEA_MASK], u[l & ALEA_MASK]);
        exp(m[(l+2) & ALEA_MASK], m[(l+1) & ALEA_MASK], u[l & ALEA_MASK]);
    }
    tim.stop();

    gmp.clear(); gmp.start();
    for (unsigned int l = 0; l < nbloops; l++) {
        mpz_powm(b[l & ALEA_MASK].get_mpz_t(),b[(l+2) & ALEA_MASK].get_mpz_t(),c[l & ALEA_MASK].get_mpz_t(), gmod.get_mpz_t());
        mpz_powm(b[(l+2) & ALEA_MASK].get_mpz_t(),b[(l+1) & ALEA_MASK].get_mpz_t(),c[l & ALEA_MASK].get_mpz_t(), gmod.get_mpz_t());
    }
    gmp.stop();


    rand(module);

    // -----------
    // Standard output for benchmark - Alexis Breust 2014/12/11
    std::cout
    << "Time: " << tim.usertime()
    << " Mflops: " << std::scientific << (double(2*nbloops))/tim.usertime()/1000.0/1000.0 << ' ' << (double(2*nbloops))/gmp.usertime()/1000.0/1000.0
    << " SIZE: " << STD_RECINT_SIZE
    << " GMP time: " << gmp.usertime()
    << ' ' << m[(int)(module )& ALEA_MASK] << ' ' << b[(int)(module)& ALEA_MASK] << std::endl ;


    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
