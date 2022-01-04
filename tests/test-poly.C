#include <iostream>
#include <givaro/givrandom.h>
#include <givaro/givtimer.h>
#include <givaro/gfq.h>
#include <givaro/givpoly1.h>
#include <givaro/givpower.h>
#include <givaro/givquotientdomain.h>

using namespace Givaro;

int main(int argc, char ** argv) {
    // argv[1] : modulo
    // argv[2] : deg max
    // argv[3] : exponent max
    // argv[4] : seed

    GFqDom<int64_t>::Residu_t MOD = (argc>1 ? (GFqDom<int64_t>::Residu_t) atoi(argv[1]) : 101);
    size_t degmax = (argc>2 ? (size_t)atoi(argv[2]) : 20);
    size_t expomax = (argc>3 ? (size_t)atoi(argv[3]) : 15);

    size_t seed = (argc>4?(size_t)atoi(argv[4]):(size_t)BaseTimer::seed());
#ifdef __GIVARO_DEBUG
    std::cerr << "seed: " << seed << std::endl;
#endif
    GivRandom generator(seed);

    GFqDom<int64_t> F(MOD);
    typedef Poly1Dom< GFqDom<int64_t>, Dense > Polynomials;
    Polynomials PD(F, Indeter("X"));

    Polynomials::Element Q;
    PD.init(Q,Degree(7));
    F.init(Q[0],7);
    F.init(Q[1],2); // Q is X^7+2X+7

    typedef QuotientDom<Polynomials> QuotRing;
    QuotRing QD(PD, Q);

    bool success = true;

        // Testing exponentiation in Quotient Domain
    uint64_t n = generator() % expomax;
    QuotRing::Element R1, R2, P;
    PD.random(generator,P,Degree(degmax));
    dom_power(R1, P, n, QD);
    PD.powmod(R2, P, Integer(n), Q);
    success &= PD.areEqual(R1,R2);
    if (! success) {
        std::cerr << "Error: " << seed << std::endl;
        PD.write(PD.write(std::cerr << "R1: ", R1) << " != ", R2) << std::endl;
    }


        // Testing divisibility in Polynomial Domain
    PD.mul(R1, P, Q);
    success &= PD.isDivisor(R1, P) && PD.isDivisor(R1, Q);
    if (! success) {
        std::cerr << "Error: " << seed << std::endl;
        PD.write(PD.write(PD.write(
            std::cerr, P) << "\n or ", Q) << "\n does not divide ", R1) << std::endl;
    }

    return (! success);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
