/*  ruint_arazi.cpp - Arazi & Qi for RecInt library test file

    Return value.
    0    No error
    != 0 Bad result for an operation

    The following constants have to be defined.
    STD_RECINT_SIZE     size of recint (> 6)
    LOOPS               number of loops
    */

#include <recint/ruint.h>
#include <recint/rmintmg.h> /* arazi_qi() */
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

int main(int argc, char ** argv)
{
    size_t nbloops = static_cast<size_t>((argc > 1)? atoi(argv[1]) : LOOPS);
    Givaro::Timer tim, gmp;

    ruint<STD_RECINT_SIZE> a[ALEA_MAX];
    ruint<STD_RECINT_SIZE> pinv[ALEA_MAX];
    mpz_class b[ALEA_MAX],c[ALEA_MAX], gmod(1);
    gmod <<= (1<<STD_RECINT_SIZE);
    //     std::cerr << "gmod: " << gmod << std::endl;

    // Randomness
    for (unsigned int i = 0; i < ALEA_MAX; i++) {
        rand(a[i]);
        if (a[i] % 2 == 0) ++a[i];
        ruint_to_mpz(b[i],a[i]);
    }

    // Random
    RecInt::srand(static_cast<unsigned long>(time(NULL)));

    tim.clear(); tim.start();
    for (UDItype l = 0; l < nbloops; l++) {
        arazi_qi(pinv[l&ALEA_MASK], a[l&ALEA_MASK]);
    }
    tim.stop();

    gmp.clear(); gmp.start();
    for (UDItype l = 0; l < nbloops; l++) {
        mpz_invert(c[l & ALEA_MASK].get_mpz_t(),b[l & ALEA_MASK].get_mpz_t(), gmod.get_mpz_t());
    }
    gmp.stop();

    ruint<STD_RECINT_SIZE> module; rand(module);

    // -----------
    // Standard output for benchmark - Alexis Breust 2014/12/11
    std::cout
    << "Time: " << tim.usertime()
    << " Mflops: " << std::scientific << (double(nbloops))/tim.usertime()/1000.0/1000.0 << ' ' << (double(nbloops))/gmp.usertime()/1000.0/1000.0
    << " SIZE: " << STD_RECINT_SIZE
    << " GMP-time: " << gmp.usertime()
    << ' ' << a[(int)(module)& ALEA_MASK] << ' ' << b[(int)(module)& ALEA_MASK] << std::endl ;

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
