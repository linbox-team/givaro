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
    size_t sP = P.size(), sQ = Q.size(), sR = R.size();
    size_t m = sP-sQ+1;
    POLY PQ;
    PD.init(PQ, Degree(sP+sQ-2));
    PD.stdmul(PQ, P, Q);
    Modular<int32_t>::Element r, pq;
    if (sR > m) return false;
    size_t i = 0;
    for (; i < sR; i++)
        if (PD.getEntry(r,i,R) != PD.getEntry(pq,i+sQ-1,PQ))
            return false;
    for (; i < m; i++)
        if (PD.getEntry(pq,i+sQ-1,PQ) != 0)
            return false;
    return true;
}

int main(int argc, char ** argv) {
    // argv[1] : modulo
    // argv[2] : value of s1 (smaller)
    // argv[3] : value of s2 (larger)
    // argv[4] : seed

    Modular<int32_t>::Residu_t MOD = (argc>1 ? (Modular<int32_t>::Residu_t) atoi(argv[1]) : 101);
    size_t s1 = (argc>2 ? (size_t)atoi(argv[2]) : 100);
    size_t s2 = (argc>3 ? (size_t)atoi(argv[3]) : s1*3/2);
    size_t seed = (argc>4?(size_t)atoi(argv[4]):(size_t)BaseTimer::seed());

#ifdef __GIVARO_DEBUG
    std::cerr << "./test-midmul " << MOD << " " << s1 << " " << s2 << " " << seed << std::endl;
#endif

    GivRandom generator(seed);

    Modular<int32_t> F(MOD);
    POLYS PD(F, Indeter("X"));

    size_t m, n;
    POLY P, Q, R;
    bool success, pass =true;

    // Balanced case

    m = s1; n = s1;
    PD.random(generator,P,Degree(m+n-2));
    PD.random(generator,Q,Degree(n-1));

        // Standard middle product

        PD.stdmidmul(R, P, Q);
        pass &= (success = testMP(PD, P, Q, R));

        if (! success) {
            std::cerr << GIV_ERROR_MSG << seed << ':' << m << ',' << n << std::endl;
            PD.write(PD.write(PD.write(
                std::cerr << "StdBalMP(P,Q) = R with" << std::endl
                          << "P = ", P) << std::endl
                          << "Q = ", Q) << std::endl
                          << "R = ", R) << std::endl;
        } else
            std::clog << "[BalStdMidMul] : " << GIV_PASSED_MSG << std::endl;

        // Karatsuba balanced middle product

        PD.karamidmul(R, P, Q);
        pass &= (success = testMP(PD, P, Q, R));

        if (! success) {
            std::cerr << GIV_ERROR_MSG << seed << ':' << m << ',' << n << std::endl;
            PD.write(PD.write(PD.write(
                std::cerr << "KaraMP(P,Q) = R with" << std::endl
                          << "P = ", P) << std::endl
                          << "Q = ", Q) << std::endl
                          << "R = ", R) << std::endl;
        } else
            std::clog << "[BalKaraMidMul] : " << GIV_PASSED_MSG << std::endl;

        // Dynamic dispatch

        PD.midmul(R, P, Q);
        pass &= (success = testMP(PD, P, Q, R));

        if (! success) {
            std::cerr << GIV_ERROR_MSG << seed << ':' << m << ',' << n << std::endl;
            PD.write(PD.write(PD.write(
                std::cerr << "BalMP(P,Q) = R with" << std::endl
                          << "P = ", P) << std::endl
                          << "Q = ", Q) << std::endl
                          << "R = ", R) << std::endl;
        } else
            std::clog << "[BalMidMul] : " << GIV_PASSED_MSG << std::endl;

    // Unbalanced case 1: m < n

    m = s1; n = s2;
    PD.random(generator,P,Degree(m+n-2));
    PD.random(generator,Q,Degree(n-1));

        // Standard middle product

        PD.stdmidmul(R, P, Q);
        pass &= (success = testMP(PD, P, Q, R));

        if (! success) {
            std::cerr << GIV_ERROR_MSG << seed << ':' << m << ',' << n << std::endl;
            PD.write(PD.write(PD.write(
                std::cerr << "UnbalStdMP1(P,Q) = R with" << std::endl
                          << "P = ", P) << std::endl
                          << "Q = ", Q) << std::endl
                          << "R = ", R) << std::endl;
        } else
            std::clog << "[UnbalStdMidMul1] : " << GIV_PASSED_MSG << std::endl;

        // Dynamic dispatch

        PD.midmul(R, P, Q);
        pass &= (success = testMP(PD, P, Q, R));

        if (! success) {
            std::cerr << GIV_ERROR_MSG << seed << ':' << m << ',' << n << std::endl;
            PD.write(PD.write(PD.write(
                std::cerr << "UnbalMP1(P,Q) = R with" << std::endl
                          << "P = ", P) << std::endl
                          << "Q = ", Q) << std::endl
                          << "R = ", R) << std::endl;
        } else
            std::clog << "[UnbalMidMul1] : " << GIV_PASSED_MSG << std::endl;

    // Unbalanced case 2: m > n

    m = s2; n = s1;
    PD.random(generator,P,Degree(m+n-2));
    PD.random(generator,Q,Degree(n-1));

        // Standard middle product

        PD.stdmidmul(R, P, Q);
        pass &= (success = testMP(PD, P, Q, R));

        if (! success) {
            std::cerr << GIV_ERROR_MSG << seed << ':' << m << ',' << n << std::endl;
            PD.write(PD.write(PD.write(
                std::cerr << "UnbalStdMP2(P,Q) = R with" << std::endl
                          << "P = ", P) << std::endl
                          << "Q = ", Q) << std::endl
                          << "R = ", R) << std::endl;
        } else
            std::clog << "[UnbalStdMidMul2] : " << GIV_PASSED_MSG << std::endl;

        // Dynamic dispatch

        PD.midmul(R, P, Q);
        pass &= (success = testMP(PD, P, Q, R));

        if (! success) {
            std::cerr << GIV_ERROR_MSG << seed << ':' << m << ',' << n << std::endl;
            PD.write(PD.write(PD.write(
                std::cerr << "UnbalMP2(P,Q) = R with" << std::endl
                          << "P = ", P) << std::endl
                          << "Q = ", Q) << std::endl
                          << "R = ", R) << std::endl;
        } else
            std::clog << "[UnbalMidMul2] : " << GIV_PASSED_MSG << std::endl;

    return (! pass);
}

#undef POLY
#undef POLYS

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
