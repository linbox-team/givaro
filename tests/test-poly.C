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

int main(int argc, char ** argv) {
    // argv[1] : modulo
    // argv[2] : deg max
    // argv[3] : exponent max
    // argv[4] : seed

    Modular<int32_t>::Residu_t MOD = (argc>1 ? (Modular<int32_t>::Residu_t) atoi(argv[1]) : 101);
    size_t degmax = (argc>2 ? (size_t)atoi(argv[2]) : 211);
    size_t expomax = (argc>3 ? (size_t)atoi(argv[3]) : 15);

    size_t seed = (argc>4?(size_t)atoi(argv[4]):(size_t)BaseTimer::seed());
#ifdef __GIVARO_DEBUG
    std::cerr << "seed: " << seed << std::endl;
#endif
    GivRandom generator(seed);

    Modular<int32_t> F(MOD);
    typedef Poly1Dom< Modular<int32_t>, Dense > Polynomials;
    Polynomials PD(F, Indeter("X"));

    Polynomials::Element Q;
    PD.init(Q,Degree(7));
    F.init(Q[0],7);
    if (MOD != 2)
        F.init(Q[1],2); // Q is X^7+2X+7
    else
        F.init(Q[1],1); // Q is X^7+X+1

    typedef QuotientDom<Polynomials> QuotRing;
    QuotRing QD(PD, Q);

    bool success, pass =true;

        // Testing exponentiation in Quotient Domain
    uint64_t n = generator() % expomax;
    QuotRing::Element R1, R2, P;
    PD.random(generator,P,Degree(degmax));
    dom_power(R1, P, n, QD);
    PD.powmod(R2, P, Integer(n), Q);
    pass &= (success = PD.areEqual(R1,R2));
    if (! success) {
        std::cerr << GIV_ERROR_MSG << seed << std::endl;
        PD.write(PD.write(std::cerr << "R1: ", R1) << " != ", R2) << std::endl;
    } else
        std::clog << "[Expo  ] : " << GIV_PASSED_MSG << std::endl;


        // Testing divisibility in Polynomial Domain
    PD.mul(R1, P, Q);
    pass &= (success = PD.isDivisor(R1, P) && PD.isDivisor(R1, Q));
    if (! success) {
        std::cerr << GIV_ERROR_MSG << seed << std::endl;
        PD.write(PD.write(PD.write(
            std::cerr, P) << "\n or ", Q) << "\n does not divide ", R1) << std::endl;
    } else
        std::clog << "[is Div] : " << GIV_PASSED_MSG << std::endl;

        // Testing modular inverses
    PD.invmod(R1, P, Q);
    PD.modin( PD.mul(R2, R1, P), Q);
    pass &= (success = PD.isOne(R2));
    if (! success) {
        std::cerr << GIV_ERROR_MSG << seed << std::endl;
        PD.write(PD.write(PD.write(PD.write(
            std::cerr << '(', R1) << ") * (", P) << ") is ", R2)
                 << " mod ", Q) << std::endl;
    } else
        std::clog << "[MdInv ] : " << GIV_PASSED_MSG << std::endl;


    PD.invmodunit(R1, P, Q);
    PD.modin( PD.mul(R2, R1, P), Q);
    pass &= (success = (PD.degree(R2) <= 0));
    if (! success) {
        std::cerr << GIV_ERROR_MSG << seed << std::endl;
        PD.write(PD.write(PD.write(PD.write(
            std::cerr << '(', R1) << ") * (", P) << ") is ", R2)
                 << " mod ", Q) << std::endl;
    } else
        std::clog << "[UMdInv] : " << GIV_PASSED_MSG << std::endl;


    PD.invmodpowx(R1, P, Degree(degmax) );
    PD.modpowxin( PD.mul(R2, R1, P), Degree(degmax));
    pass &= (success = (PD.degree(R2) <= 0));
    if (! success) {
        std::cerr << GIV_ERROR_MSG << seed << ':' << PD.degree(R2) << std::endl;
        PD.write(PD.write(PD.write(
            std::cerr << "Rem((", R1) << ") * (", P) << "),X^" << degmax << ",X) mod " << MOD << " = ", R2) << ';' << std::endl;
    } else
        std::clog << "[MdInvX] : " << GIV_PASSED_MSG << std::endl;



    return (! pass);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
