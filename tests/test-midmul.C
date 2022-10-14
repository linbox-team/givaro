#include <iostream>
#include <givaro/givrandom.h>
#include <givaro/givtimer.h>
#include <givaro/modular.h>
#include <givaro/givpoly1.h>
#include <givaro/givpower.h>
#include <givaro/givquotientdomain.h>

using namespace Givaro;

#define GIV_PASSED_MSG "\033[1;32mPASSED.\033[0m"
#define GIV_ERROR_MSG "\033[1;31m****** ERROR ******\033[0m "

#define POLYS Poly1Dom< Modular<int32_t>, Dense >
#define POLY POLYS::Element

bool testMP(POLYS PD, POLY P, POLY Q, POLY R) {
    POLY PQ;
    PD.stdmul(PQ, P, Q);
    size_t sP = P.size(), sQ = Q.size(), sR = R.size();
    size_t m = sP-sQ+1;
    Modular<int32_t>::Element r, pq;
    if (sR != m) return false;
    for (size_t i = 0; i < sR; i++)
        if (PD.getEntry(r,i,R) != PD.getEntry(pq,i+sQ-1,PQ))
            return false;
    return true;
}

int main(int argc, char ** argv) {
    // argv[1] : modulo
    // argv[2] : value of m
    // argv[3] : value of n
    // argv[4] : seed

    Modular<int32_t>::Residu_t MOD = (argc>1 ? (Modular<int32_t>::Residu_t) atoi(argv[1]) : 101);
    size_t m = (argc>2 ? (size_t)atoi(argv[2]) : 100);
    size_t n = (argc>3 ? (size_t)atoi(argv[3]) : 100);

    size_t seed = (argc>4?(size_t)atoi(argv[4]):(size_t)BaseTimer::seed());
#ifdef __GIVARO_DEBUG
    std::cerr << "seed: " << seed << std::endl;
#endif
    GivRandom generator(seed);

    Modular<int32_t> F(MOD);
    POLYS PD(F, Indeter("X"));

    POLY P, Q, R;
    PD.random(generator,P,Degree(m+n));
    PD.random(generator,Q,Degree(n-1));
    bool success, pass =true;



    // Standard middle product

    PD.stdmidmul(R, P, Q);
    pass &= (success = testMP(PD, P, Q, R));

    if (! success) {
        std::cerr << GIV_ERROR_MSG << seed << ':' << m << ',' << n << std::endl;
        PD.write(PD.write(PD.write(
            std::cerr << "MP((", P) << "), (", Q) << ") = ", R) << ");" << std::endl;
    } else
        std::clog << "[stdMidMul] : " << GIV_PASSED_MSG << std::endl;



    return (! pass);
}

#undef POLY
#undef POLYS

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
