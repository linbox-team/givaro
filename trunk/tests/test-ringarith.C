// Copyright(c)'1994-2009,2010 by The Givaro group
// This file is part of Givaro.
// BB: adapted,enhanced from examples/FiniteFields/ff-arith.C
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

#include <iostream>
#include <givaro/modular.h>
#include <givaro/montgomery.h>
#include <givaro/givpoly1.h>
#include <givaro/givinteger.h>

using namespace Givaro;

#define TESTE_EG( a, b )						\
    if (!F.areEqual((a),(b))) {						\
	F.write( F.write(std::cout,a) << "!=",b)			\
	    << " failed (at line " <<  __LINE__ << ")" << std::endl;	\
	return -1;							\
    }

#define JETESTE( a, s )					\
    if (TestRing( (a), (s)) ) {				\
	std::cout << #a << " failed !" << std::endl;	\
	return -1;					\
    }

#define JEPOLTESTE( a, s )				\
    if (TestPolRing( (a), (s) ) ) {			\
	std::cout << #a << " failed !" << std::endl;	\
	return -1;					\
    }

#define JEONETESTE( F, x, y )						\
    if (TestOneRing(F,x,y)) {						\
	std::cout << #x << " " << #y << " failed !" << std::endl;	\
	return -1;							\
    }

template<class Ring>
int TestOneRing(const Ring& F, const typename Ring::Element& x, const typename Ring::Element& y)
{
#ifdef GIVARO_DEBUG
    std::cerr << "Testing " ;
    F.write(std::cerr) << " : " << std::endl;
#endif

    typename Ring::Element a, b, c, d,a_,b_,c_,d_;
    typename Ring::Element e,e_;

    F.init(a, 0UL);
    TESTE_EG(a, F.zero);
    F.init(a, 1UL);
    //         F.write(std::cerr) << std::endl;
    //         F.write(std::cerr << "a: ", a) << std::endl;
    //         F.write(std::cerr << "1: ", F.one) << std::endl;
    TESTE_EG(a, F.one);

    F.assign(a, x);
    F.assign(b, y);
    F.init(c);            // empty constructor
    F.init(d);            // empty constructor
	
    F.add(c, a, b);       // c = a+b
    F.assign(c_,c);       // c_ <- c
    TESTE_EG(c,c_);
	
    F.subin(c_,a);
    //         F.write(std::cerr) << std::endl;
    //         F.write(std::cerr << "a: ", a) << std::endl;
    //         F.write(std::cerr << "b: ", b) << std::endl;
    //         F.write(std::cerr << "c: ", c) << std::endl;
    //         F.write(std::cerr << "c_: ", c_) << std::endl;
    TESTE_EG(b,c_);

    F.axpy(d, a, b, c); // d = a*b + c;
    F.init(d_);
    F.axmy(d_,a,b,c); // d_ = a*b - c
    F.addin(d_,c);
    F.subin(d,c);
    //         F.write(std::cerr) << std::endl;
    //         F.write(std::cerr << "a: ", a) << std::endl;
    //         F.write(std::cerr << "b: ", b) << std::endl;
    //         F.write(std::cerr << "c: ", c) << std::endl;
    //         F.write(std::cerr << "d: ", d) << std::endl;
    //         F.write(std::cerr << "d_: ", d_) << std::endl;
    TESTE_EG(d_,d);

    F.axpy(d, a, b, c); // d = a*b + c;
    F.init(d_);
    F.assign(d_,c);
    F.axpyin(d_,a,b);
    TESTE_EG(d_, d);

    F.sub(d,a,b); // d = a -b
    F.add(c,a,b); // c = a+b
    F.init(e);
    F.init(e_);
    F.mul(e,d,c); // e = d*c;
    F.mul(a_,a,a); // a_ = a*a
    F.mul(b_,b,b); // b_ = b*b
    F.sub(e_,a_,b_); // e_ = a_ - b_
    TESTE_EG(e,e_); // a^2 - b^2 = (a-b)(a+b)

    F.maxpy(e, a, b, d); // e = d-a*b
    F.assign(e_,d);
    F.maxpyin(e_, a, b); // e = d - a*b;
    TESTE_EG(e,e_);

    F.axmy(e, a, b, d); // e = a*b -d;
    F.assign(e_,d);
    F.maxpyin(e_, a, b); // e = d - a*b;
    F.negin(e_);
    TESTE_EG(e,e_);

    F.maxpy(e, a, b, d); // e = d-a*b;
    F.assign(e_,d);
    F.axmyin(e_, a, b); // e_ = a*b-e_ = a*b-d
    F.negin(e_);
    //         F.write(std::cerr) << std::endl;
    //         F.write(std::cerr << "a: ", a) << std::endl;
    //         F.write(std::cerr << "b: ", b) << std::endl;
    //         F.write(std::cerr << "c: ", c) << std::endl;
    //         F.write(std::cerr << "d: ", d) << std::endl;
    //         F.write(std::cerr << "e: ", e) << std::endl;
    //         F.write(std::cerr << "e_: ", e_) << std::endl;
    TESTE_EG(e,e_);

#ifdef GIVARO_DEBUG
    F.write(std::cerr );
    std::cerr  << " done." << std::endl;
#endif
	
    return 0;
}

#ifndef NBITERD
#define NBITER 50
#endif

template<class Ring>
int TestRing(const Ring& F, const unsigned  long seed)
{
    typename Ring::Element x, y;
    typename Ring::RandIter g(F, seed);
    
    F.init(x, 7UL);
    F.init(y, -29.3);
    JEONETESTE(F,x,y);
    
    for (size_t i = 0; i< NBITER; ++i) {
	g.random(x); g.random(y);
        JEONETESTE(F,x,y);
    }
    
    return 0;
}

#ifndef DEGMAX
#define DEGMAX 75
#endif
#ifndef NBITERD
#define NBITERD 10
#endif

template<class Ring>
int TestPolRing(const Ring& F, const unsigned long seed)
{
    GivRandom generator(seed);
    srand48((long)seed);

    for (size_t i = 0; i < NBITERD; ++i) {
        int d1 = int (lrand48() % DEGMAX);
        int d2 = int (lrand48() % DEGMAX);
        typename Ring::Element x, d, z, o;
	
        do { F.random(generator, x, Degree(d1)); } while(F.isZero(x));
        do { F.random(generator, d, Degree(d2)); } while(F.isZero(d));
        JEONETESTE(F,x,d);
	
        do { F.random(generator, z, Degree(0)); } while(F.isZero(z));
        JEONETESTE(F,x,z);
        JEONETESTE(F,z,x);
	
        do { F.random(generator, o, Degree(1)); } while(F.isZero(o));
        JEONETESTE(F,d,o);
        JEONETESTE(F,o,d);
    }
    
    return 0;
}

int main(int argc, char ** argv)
{
    auto seed = static_cast<unsigned long>(argc>1?atoi(argv[1]):BaseTimer::seed());
#ifdef GIVARO_DEBUG
    std::cerr << "seed: " << seed << std::endl;
#endif
    Integer::seeding(seed);
    RecInt::srand(seed);

    //-------------//
    //----- 2 -----//

    // modulo 2 over 8 bits
    Modular<int8_t> C2(2); JETESTE(C2,seed);
    Modular<uint8_t> UC2(2); JETESTE(UC2,seed);

    // modulo 2 over 16 bits
    Modular<int16_t> S2(2); JETESTE(S2,seed);
    Modular<uint16_t> US2(2); JETESTE(US2,seed);

    // modulo 2 over 32 bits
    Modular<int32_t> Z2(2); JETESTE(Z2,seed);
    Modular<uint32_t> UZ2(2); JETESTE(UZ2,seed);

    // modulo 2 over 64 bits
    Modular<int64_t> LL2(2); JETESTE(LL2,seed);
    Modular<uint64_t> ULL2(2); JETESTE(ULL2,seed);

    // modulo 2 fully tabulated
    Modular<Log16> L2(2); JETESTE(L2,seed);

    // modulo 2 over arbitrary size
    Modular<Integer> IntZ2(2); JETESTE(IntZ2,seed);

    // modulo 2 over 128 bits
    Modular<RecInt::ruint128> UR2(2); JETESTE(UR2,seed);
    // Not odd Montgomery<RecInt::ruint128> UM2(2); JETESTE(UM2,seed);

    //-------------//
    //----- 4 -----//

    // modulo 4 over 8 bits
    Modular<int8_t> C4(4); JETESTE(C4,seed);
    Modular<uint8_t> UC4(4); JETESTE(UC4,seed);

    // modulo 4 over 16 bits
    Modular<int16_t> S4(4); JETESTE(S4,seed);
    Modular<uint16_t> US4(4); JETESTE(US4,seed);

    // modulo 4 over 32 bits
    Modular<int32_t> Z4(4); JETESTE(Z4,seed);
    Modular<uint32_t> UZ4(4); JETESTE(UZ4,seed);

    // modulo 4 over 64 bits
    Modular<int64_t> LL4(4); JETESTE(LL4,seed);
    Modular<uint64_t> ULL4(4); JETESTE(ULL4,seed);

    // modulo 4 fully tabulated
    // Not prime Modular<Log16> L4(4); JETESTE(L4,seed);

    // modulo 4 over arbitrary size
    Modular<Integer> IntZ4(4); JETESTE(IntZ4,seed);

    // modulo 4 over 128 bits
    Modular<RecInt::ruint128> UR4(4); JETESTE(UR4,seed);
    // Not odd Montgomery<RecInt::ruint128> UM4(4); JETESTE(UM4,seed);

    //--------------//
    //----- 13 -----//

    // modulo 13 over 8 bits
    Modular<int8_t> C13(13); JETESTE(C13,seed);
    Modular<uint8_t> UC13(13); JETESTE(UC13,seed);

    // modulo 13 over 16 bits
    Modular<int16_t> S13(13); JETESTE(S13,seed);
    Modular<uint16_t> US13(13); JETESTE(US13,seed);

    // modulo 13 over 32 bits
    Modular<int32_t> Z13(13); JETESTE(Z13,seed);
    Modular<uint32_t> UZ13(13); JETESTE(UZ13,seed);

    // modulo 13 over 64 bits
    Modular<int64_t> LL13(13); JETESTE(LL13,seed);
    Modular<uint64_t> ULL13(13); JETESTE(ULL13,seed);

    // modulo 13 fully tabulated
    Modular<Log16> L13(13); JETESTE(L13,seed);

    // modulo 13 over arbitrary size
    Modular<Integer> IntZ13(13); JETESTE(IntZ13,seed);

    // modulo 13 over 128 bits
    Modular<RecInt::ruint128> UR13(13); JETESTE(UR13,seed);
    Montgomery<RecInt::ruint128> UM13(13); JETESTE(UM13,seed);

    // Polynomial
    Poly1Dom< Modular<int8_t>, Dense > CP13(C13, "X"); JEPOLTESTE(CP13,seed);
    Poly1Dom< Modular<uint8_t>, Dense > UCP13(UC13, "X"); JEPOLTESTE(UCP13,seed);
    Poly1Dom< Modular<int16_t>, Dense > SP13(S13, "X"); JEPOLTESTE(SP13,seed);
    Poly1Dom< Modular<uint16_t>, Dense > USP13(US13, "X"); JEPOLTESTE(USP13,seed);
    Poly1Dom< Modular<int32_t>, Dense > ZP13(Z13, "X"); JEPOLTESTE(ZP13,seed);
    Poly1Dom< Modular<uint32_t>, Dense > UZP13(UZ13, "X"); JEPOLTESTE(UZP13,seed);
    Poly1Dom< Modular<int64_t>, Dense > LLP13(LL13, "X"); JEPOLTESTE(LLP13,seed);
    Poly1Dom< Modular<uint64_t>, Dense > ULLP13(ULL13, "X"); JEPOLTESTE(ULLP13,seed);
    Poly1Dom< Modular<Integer>, Dense > IntZP13(IntZ13, "X"); JEPOLTESTE(IntZP13,seed);
    Poly1Dom< Modular<RecInt::ruint128>, Dense > URP13(UR13, "X"); JEPOLTESTE(URP13,seed);

    //--------------//
    //----- 75 -----//

    // modulo 75 over 8 bits
    Modular<int8_t> C75(75); JETESTE(C75,seed);
    Modular<uint8_t> UC75(75); JETESTE(UC75,seed);

    // modulo 75 over 16 bits
    Modular<int16_t> S75(75); JETESTE(S75,seed);
    Modular<uint16_t> US75(75); JETESTE(US75,seed);

    // modulo 75 over 32 bits
    Modular<int32_t> Z75(75); JETESTE(Z75,seed);
    Modular<uint32_t> UZ75(75); JETESTE(UZ75,seed);

    // modulo 75 over 64 bits
    Modular<int64_t> LL75(75); JETESTE(LL75,seed);
    Modular<uint64_t> ULL75(75); JETESTE(ULL75,seed);

    // modulo 75 fully tabulated
    // Not prime Modular<Log16> L75(75); JETESTE(L75,seed);

    // modulo 75 over arbitrary size
    Modular<Integer> IntZ75(75); JETESTE(IntZ75,seed);

    // modulo 75 over 128 bits
    Modular<RecInt::ruint128> UR75(75); JETESTE(UR75,seed);
    Montgomery<RecInt::ruint128> UM75(75); JETESTE(UM75,seed);

    // Polynomial
    Poly1Dom< Modular<int8_t>, Dense > CP75(C75, "X"); JEPOLTESTE(CP75,seed);
    Poly1Dom< Modular<uint8_t>, Dense > UCP75(UC75, "X"); JEPOLTESTE(UCP75,seed);
    Poly1Dom< Modular<int16_t>, Dense > SP75(S75, "X"); JEPOLTESTE(SP75,seed);
    Poly1Dom< Modular<uint16_t>, Dense > USP75(US75, "X"); JEPOLTESTE(USP75,seed);
    Poly1Dom< Modular<int32_t>, Dense > ZP75(Z75, "X"); JEPOLTESTE(ZP75,seed);
    Poly1Dom< Modular<uint32_t>, Dense > UZP75(UZ75, "X"); JEPOLTESTE(UZP75,seed);
    Poly1Dom< Modular<int64_t>, Dense > LLP75(LL75, "X"); JEPOLTESTE(LLP75,seed);
    Poly1Dom< Modular<uint64_t>, Dense > ULLP75(ULL75, "X"); JEPOLTESTE(ULLP75,seed);
    Poly1Dom< Modular<Integer>, Dense > IntZP75(IntZ75, "X"); JEPOLTESTE(IntZP75,seed);
    Poly1Dom< Modular<RecInt::ruint128>, Dense > URP75(UR75, "X"); JEPOLTESTE(URP75,seed);

    // Inception
    Poly1Dom< Poly1Dom< Modular<Integer>, Dense >, Dense> IntZPP75(IntZP75, "Y");
    JEPOLTESTE(IntZPP75,seed);

    Poly1Dom< Poly1Dom< Poly1Dom< Modular<Integer>, Dense >, Dense>, Dense > IntZPPP75(IntZPP75, "Z");
    JEPOLTESTE(IntZPPP75,seed);

    return 0;
}
