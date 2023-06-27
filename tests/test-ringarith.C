// Copyright(c)'1994-2009,2010 by The Givaro group
// This file is part of Givaro.
// BB: adapted,enhanced from examples/FiniteFields/ff-arith.C
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

#include <iostream>
#include <givaro/modular.h>
#include <givaro/modular-balanced.h>
#include <givaro/modular-extended.h>
//#include <givaro/montgomery.h>
#include <givaro/givpoly1.h>
#include <givaro/givinteger.h>
#include <givaro/zring.h>
#include <givaro/gfq.h>
#include <recint/recint.h>
//#include <givaro/modular-extended.h>
//#include <givaro/modular-general.h>
//#include <givaro/modular-integral.h>
//#include <givaro/modular-floating.h>
//#include <givaro/modular-integer.h>
//#include <givaro/modular-ruint.h>
//#include <givaro/modular-log16.h>
//#include <givaro/modular-inttype.h>

using namespace Givaro;


#define TESTE_EG( a, b )						\
if (!F.areEqual((a),(b))) {						\
    F.write( F.write(std::cout,a) << "!=",b)			\
    << " failed (at line " <<  __LINE__ << ")" << std::endl;	\
    return -1;							\
}

#define TESTE_T( b )						\
if (!b) {						\
    F.write(std::cout)		\
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
#ifdef __GIVARO_DEBUG
    std::cerr << "Testing " ;
    F.write(std::cerr) << " : " << std::endl;
#endif

    typename Ring::Element a, b, c, d,a_,b_,c_,d_;
    typename Ring::Element e,e_;

    F.init(a, 0U);
    TESTE_EG(a, F.zero);
    TESTE_T(F.isZero(a));

    F.init(a, 1U);
    //         F.write(std::cerr) << std::endl;
    //         F.write(std::cerr << "a: ", a) << std::endl;
    //         F.write(std::cerr << "1: ", F.one) << std::endl;
    TESTE_EG(a, F.one);
    TESTE_T(F.isOne(a));
    TESTE_T(F.isUnit(a));

    F.init(a, -1);
    TESTE_EG(a, F.mOne);
    TESTE_T(F.isMOne(a));

    F.assign(a, x);
    F.assign(b, y);
    F.init(c);            // empty constructor
    F.init(d);            // empty constructor

    F.add(c, a, b);       // c = a+b
    F.assign(c_,c);       // c_ <- c
    TESTE_EG(c,c_);

    F.subin(c_,a);
    //F.write(std::cerr) << std::endl;
    //F.write(std::cerr << "a: ", a) << std::endl;
    //F.write(std::cerr << "b: ", b) << std::endl;
    //F.write(std::cerr << "c: ", c) << std::endl;
    //F.write(std::cerr << "c_: ", c_) << std::endl;
    TESTE_EG(b,c_);

    F.axpy(d, a, b, c); // d = a*b + c;
    F.init(d_);
    F.axmy(d_,a,b,c); // d_ = a*b - c
    F.addin(d_,c);
    F.subin(d,c);
    //F.write(std::cerr) << std::endl;
    //F.write(std::cerr << "a: ", a) << std::endl;
    //F.write(std::cerr << "b: ", b) << std::endl;
    //F.write(std::cerr << "c: ", c) << std::endl;
    //F.write(std::cerr << "d: ", d) << std::endl;
    //F.write(std::cerr << "d_: ", d_) << std::endl;
    TESTE_EG(d_,d);

    F.axpy(d, a, b, c); // d = a*b + c;
    F.init(d_);
    F.assign(d_,c);
    F.axpyin(d_,a,b);
    //F.write(std::cerr) << std::endl;
    //F.write(std::cerr << "a: ", a) << std::endl;
    //F.write(std::cerr << "b: ", b) << std::endl;
    //F.write(std::cerr << "c: ", c) << std::endl;
    //F.write(std::cerr << "d: ", d) << std::endl;
    //F.write(std::cerr << "d_: ", d_) << std::endl;
    TESTE_EG(d_, d);

    F.sub(d,a,b); // d = a -b
    F.add(c,a,b); // c = a+b
    F.init(e);
    F.init(e_);
    F.mul(e,d,c); // e = d*c;
    F.mul(a_,a,a); // a_ = a*a
    F.mul(b_,b,b); // b_ = b*b
    F.sub(e_,a_,b_); // e_ = a_ - b_
    //         F.write(std::cerr) << std::endl;
    //         F.write(std::cerr << "a: ", a) << std::endl;
    //         F.write(std::cerr << "a_: ", a_) << std::endl;
    //         F.write(std::cerr << "b: ", b) << std::endl;
    //         F.write(std::cerr << "b_: ", b_) << std::endl;
    //         F.write(std::cerr << "c: ", c) << std::endl;
    //         F.write(std::cerr << "d: ", d) << std::endl;
    //         F.write(std::cerr << "e: ", e) << std::endl;
    //         F.write(std::cerr << "e_: ", e_) << std::endl;
    TESTE_EG(e,e_); // a^2 - b^2 = (a-b)(a+b)

    F.maxpy(e, a, b, d); // e = d-a*b
    F.assign(e_,d);
    F.maxpyin(e_, a, b); // e = d - a*b;
    TESTE_EG(e,e_);

    F.assign(d,a);
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

    F.init(a, 3);
    F.assign(b , y);
    F.mul(c,a,b);
    F.subin(c,b);
    F.subin(c,b);
    F.subin(c,b);
    F.init(d, 0U);
    TESTE_EG(c, d);

    //F.write(std::cerr) << std::endl;
    F.init(a, -3);          // F.write(std::cerr << "a: ", a) << std::endl;
    F.assign(b , y);        // F.write(std::cerr << "b: ", b) << std::endl;
    F.mul(c,a,b);           // F.write(std::cerr << "c: ", c) << std::endl;
    F.addin(c,b);           // F.write(std::cerr << "c: ", c) << std::endl;
    F.addin(c,b);           // F.write(std::cerr << "c: ", c) << std::endl;
    F.addin(c,b);           // F.write(std::cerr << "c: ", c) << std::endl;
    F.init(d, 0U);

    //F.write(std::cerr << "d: ", d) << std::endl;


    TESTE_EG(c, d);

#ifdef __GIVARO_DEBUG
    F.write(std::cerr );
    std::cerr  << " done." << std::endl;
#endif

    return 0;
}

#ifndef NBITERD
#define NBITER 50
#endif

template<class Ring>
int TestRing(const Ring& F, const uint64_t seed)
{
    typename Ring::Element x, y;
    typename Ring::RandIter g(F, seed);

    F.init(x, 7U);
    F.init(y, -29.0);
    JEONETESTE(F,x,y);

    F.init(x, Ring::maxCardinality()-1U);
    F.init(y, Ring::maxCardinality()-1U);
    JEONETESTE(F,x,y);

    F.assign(x, F.maxElement());
    F.assign(y, F.maxElement());
    JEONETESTE(F,x,y);

    F.assign(x, F.minElement());
    F.assign(y, F.maxElement());
    JEONETESTE(F,x,y);

    F.assign(x, F.minElement());
    F.assign(y, F.minElement());
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
int TestPolRing(const Ring& F, const uint64_t seed)
{
    GivRandom generator(seed);
    srand48((int64_t)seed);

    for (size_t i = 0; i < NBITERD; ++i) {
	int64_t d1 = int64_t (lrand48() % DEGMAX);
	int64_t d2 = int64_t (lrand48() % DEGMAX);
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

template<class Ring>
int TestInv(const Ring& F, const uint64_t seed)
{
    typename Ring::Element a, inva, invinva;
    typename Ring::RandIter gg(F, seed);
    typename Ring::NonZeroRandIter g(gg);

    F.init(a);
    F.init(inva);
    F.init(invinva);

    for (int i=0; i < 100; i++)
    {
	g.random(a);
	F.inv(inva, a);
	F.inv(invinva, inva);

	// F.write(std::cerr);
	// F.write(std::cerr << " => a: ", a);
	// F.write(std::cerr << "; inva: ", inva);
	// F.write(std::cerr << "; invinva: ", invinva) << std::endl;

	TESTE_EG(a, invinva);
    }


    return 0;
}


int main(int argc, char ** argv)
{
    auto seed = static_cast<uint64_t>(argc>1?atoi(argv[1]):BaseTimer::seed());
#ifdef __GIVARO_DEBUG
    std::cerr << "seed: " << seed << std::endl;
#endif
    Integer::seeding(seed);
    RecInt::srand(seed);

    using ModularCUS = Modular<int8_t, uint16_t>;
    using ModularSUZ = Modular<int16_t, uint32_t>;
    using ModularZULL = Modular<int32_t, uint64_t>;
    using ModularUCUS = Modular<uint8_t, uint16_t>;
    using ModularUSUZ = Modular<uint16_t, uint32_t>;
    using ModularUZULL = Modular<uint32_t, uint64_t>;
    using ModularFD = Modular<float, double>;

#ifdef __GIVARO_HAVE_INT128
    using ModularLLULLL = Modular<int64_t, uint128_t>;
    using ModularULLULLL = Modular<uint64_t, uint128_t>;
#endif

#define TEST_SPECIFIC(Ring, Name, Modulus...)	\
    std::cout << "TEST_SPECIFIC: " << #Name << std::endl; \
    Ring Name(Modulus);				\
    JETESTE(Name, seed);

    //Modular<float> GF5(5);
    //float elt, inv;
    //GF5.init(elt, 2);
    //GF5.inv(inv, elt);
    //GF5.write(std::cout, inv) << " ####\n";

    //-------------//
    //----- 4 -----//

    TEST_SPECIFIC(Modular<int8_t>, C4, 4);
    TEST_SPECIFIC(Modular<int16_t>, S4, 4);
    TEST_SPECIFIC(Modular<int32_t>, Z4, 4);
    TEST_SPECIFIC(Modular<int64_t>, LL4, 4);
    TEST_SPECIFIC(Modular<uint8_t>, UC4, 4);
    TEST_SPECIFIC(Modular<uint16_t>, US4, 4);
    TEST_SPECIFIC(Modular<uint32_t>, UZ4, 4);
    TEST_SPECIFIC(Modular<uint64_t>, ULL4, 4);
    TEST_SPECIFIC(ModularCUS, CUS4, 4);
    TEST_SPECIFIC(ModularSUZ, SUZ4, 4);
    TEST_SPECIFIC(ModularZULL, ZULL4, 4);
    TEST_SPECIFIC(ModularUCUS, UCUS4, 4);
    TEST_SPECIFIC(ModularUSUZ, USUZ4, 4);
    TEST_SPECIFIC(ModularUZULL, UZULL4, 4);
    TEST_SPECIFIC(ModularFD, FD4, 4);
#ifdef __GIVARO_HAVE_INT128
    TEST_SPECIFIC(ModularLLULLL, LLULLL4, 4);
    TEST_SPECIFIC(ModularULLULLL, ULLULLL4, 4);
#endif

    TEST_SPECIFIC(Modular<float>, F4, 4);
    TEST_SPECIFIC(Modular<double>, D4, 4);
    TEST_SPECIFIC(ModularExtended<float>, MEF4, 4);
    TEST_SPECIFIC(ModularExtended<double>, MED4, 4);
    TEST_SPECIFIC(Modular<Integer>, I4, 4);
    TEST_SPECIFIC(Modular<RecInt::ruint128>, RU4, 4);
    TEST_SPECIFIC(Modular<RecInt::rint128>, R4, 4);
    TEST_SPECIFIC(ZRing<Integer>, ZR4, 4);

    TEST_SPECIFIC(Modular<Log16>, L5, 5);
    TEST_SPECIFIC(Modular<Log16>, L7, 7);

    //--------------//
    //----- 75 -----//

    TEST_SPECIFIC(Modular<Log16>, L79, 79);

    TEST_SPECIFIC(Modular<int8_t>, C75, 13);
    TEST_SPECIFIC(Modular<int16_t>, S75, 75);
    TEST_SPECIFIC(Modular<int32_t>, Z75, 75);
    TEST_SPECIFIC(Modular<int64_t>, LL75, 75);
    TEST_SPECIFIC(Modular<uint8_t>, UC75, 13);
    TEST_SPECIFIC(Modular<uint16_t>, US75, 75);
    TEST_SPECIFIC(Modular<uint32_t>, UZ75, 75);
    TEST_SPECIFIC(Modular<uint64_t>, ULL75, 75);
    TEST_SPECIFIC(ModularCUS, CUS75, 75);
    TEST_SPECIFIC(ModularSUZ, SUZ75, 75);
    TEST_SPECIFIC(ModularZULL, ZULL75, 75);
    TEST_SPECIFIC(ModularUCUS, UCUS75, 75);
    TEST_SPECIFIC(ModularUSUZ, USUZ75, 75);
    TEST_SPECIFIC(ModularUZULL, UZULL75, 75);
    TEST_SPECIFIC(ModularFD, FD75, 75);
#ifdef __GIVARO_HAVE_INT128
    TEST_SPECIFIC(ModularLLULLL, LLULLL75, 75);
    TEST_SPECIFIC(ModularULLULLL, ULLULLL75, 75);
#endif


    TEST_SPECIFIC(Modular<float>, F75, 75);
    TEST_SPECIFIC(Modular<double>, D75, 75);
    TEST_SPECIFIC(ModularExtended<float>, MEF75, 75);
    TEST_SPECIFIC(ModularExtended<double>, MED75, 75);
    TEST_SPECIFIC(Modular<Integer>, I75, 75);
    TEST_SPECIFIC(Modular<RecInt::ruint128>, RU75, 75);
    TEST_SPECIFIC(Modular<RecInt::rint128>, R75, 75);

    TEST_SPECIFIC(ModularBalanced<int32_t>, BZ75, 75);
    TEST_SPECIFIC(ModularBalanced<int64_t>, BLL75, 75);
    TEST_SPECIFIC(ModularBalanced<float>, BF75, 75);
    TEST_SPECIFIC(ModularBalanced<double>, BD75, 75);
    //TEST_SPECIFIC(Montgomery<int32_t>, MZ75, 75);
    //TEST_SPECIFIC(Montgomery<RecInt::ruint128>, MRU75, 75);

#define TEST_POLYNOMIAL(BaseRing, Name, BaseRingName)		\
    Poly1Dom<BaseRing, Dense> Name(BaseRingName, "X");		\
    JEPOLTESTE(Name, seed);

    TEST_POLYNOMIAL(Modular<Log16>, PL79, L79);
    TEST_POLYNOMIAL(Modular<int8_t>, PC75, C75);
    TEST_POLYNOMIAL(Modular<int16_t>, PS75, S75);
    TEST_POLYNOMIAL(Modular<int32_t>, PZ75, Z75);
    TEST_POLYNOMIAL(Modular<int64_t>, PLL75, LL75);
    TEST_POLYNOMIAL(Modular<uint8_t>, PUC75, UC75);
    TEST_POLYNOMIAL(Modular<uint16_t>, PUS75, US75);
    TEST_POLYNOMIAL(Modular<uint32_t>, PUZ75, UZ75);
    TEST_POLYNOMIAL(Modular<uint64_t>, PULL75, ULL75);
    TEST_POLYNOMIAL(ModularCUS, PCUS75, CUS75);
    TEST_POLYNOMIAL(ModularSUZ, PSUZ75, SUZ75);
    TEST_POLYNOMIAL(ModularZULL, PZULL75, ZULL75);
    TEST_POLYNOMIAL(ModularUCUS, PUCUS75, UCUS75);
    TEST_POLYNOMIAL(ModularUSUZ, PUSUZ75, USUZ75);
    TEST_POLYNOMIAL(ModularUZULL, PUZULL75, UZULL75);
    TEST_POLYNOMIAL(ModularFD, PFD75, FD75);
    TEST_POLYNOMIAL(Modular<float>, PF75, F75);
    TEST_POLYNOMIAL(Modular<double>, PD75, D75);
    TEST_POLYNOMIAL(ModularExtended<float>, PMEF75, MEF75);
    TEST_POLYNOMIAL(ModularExtended<double>, PMED75, MED75);
    TEST_POLYNOMIAL(Modular<Integer>, PI75, I75);
    TEST_POLYNOMIAL(Modular<RecInt::ruint128>, PRU75, RU75);
    TEST_POLYNOMIAL(Modular<RecInt::rint128>, PR75, R75);

    TEST_POLYNOMIAL(ModularBalanced<int32_t>, MBZ75, BZ75);
    TEST_POLYNOMIAL(ModularBalanced<int64_t>, MBLL75, BLL75);
    TEST_POLYNOMIAL(ModularBalanced<float>, MBF75, BF75);
    TEST_POLYNOMIAL(ModularBalanced<double>, MBD75, BD75);
    //TEST_POLYNOMIAL(Montgomery<int32_t>, PMZ75, MZ75);
    // @bug Convert to double inside? //TEST_POLYNOMIAL(Montgomery<RecInt::ruint128>, PMRU75, MRU75);

    //TEST_POLYNOMIAL(decltype(PI75), PPI75, PI75);
    //TEST_POLYNOMIAL(decltype(PPI75), PPPI75, PPI75);

    //--------------------------------//
    //----- Modulo maximal prime -----//

#define TEST_LAST(Field, Name)		\
    std::cout << "TEST_LAST: " << #Name; \
    Field Name(Field::maxCardinality());	\
    Name.write(std::cout << " (", Field::maxCardinality()) << ")"<< std::endl; \
    JETESTE(Name, seed);

    TEST_LAST(Modular<Log16>, Lmax);
    TEST_LAST(Modular<int8_t>, Cmax);
    TEST_LAST(Modular<int16_t>, Smax);
    TEST_LAST(Modular<int32_t>, Zmax);
    TEST_LAST(Modular<int64_t>, LLmax);
    TEST_LAST(Modular<uint8_t>, UCmax);
    TEST_LAST(Modular<uint16_t>, USmax);
    TEST_LAST(Modular<uint32_t>, UZmax);
    TEST_LAST(Modular<uint64_t>, ULLmax);
    TEST_LAST(ModularCUS, CUSmax);
    TEST_LAST(ModularSUZ, SUZmax);
    TEST_LAST(ModularZULL, ZULLmax);
    TEST_LAST(ModularUCUS, UCUSmax);
    TEST_LAST(ModularUSUZ, USUZmax);
    TEST_LAST(ModularUZULL, UZULLmax);
#ifdef __GIVARO_HAVE_INT128
    TEST_LAST(ModularLLULLL, LLULLLmax);
    TEST_LAST(ModularULLULLL, ULLULLLmax);
#endif

    TEST_LAST(Modular<float>, Fmax);
    TEST_LAST(Modular<double>, Dmax);
    TEST_LAST(ModularExtended<float>, MEFmax);
    TEST_LAST(ModularExtended<double>, MEDmax);
    TEST_LAST(ModularFD, FDmax);
    TEST_LAST(Modular<RecInt::ruint128>, RUmax);
    TEST_LAST(Modular<RecInt::rint128>,  Rmax);
    typedef Modular<RecInt::ruint128,RecInt::ruint256> MyMod;
    TEST_LAST(MyMod, RUExtmax);
    TEST_LAST(ModularBalanced<int32_t>, BZmax);
    TEST_LAST(ModularBalanced<int64_t>, BLLmax);
    TEST_LAST(ModularBalanced<float>, BFmax);
    TEST_LAST(ModularBalanced<double>, BDmax);

    //TEST_LAST(Montgomery<int32_t>, MZmax);
    //TEST_LAST(Montgomery<RecInt::ruint128>, MRUmax);

    // -----------------------
    //  Test inversions
    // -----------------------


#define TEST_INV(Field, Name, Prime) \
    std::cout << "TEST_INV: " << #Name; \
    Field Name(Prime); \
    std::cout << " (" << (Integer)Name.cardinality() << ',' << (Integer)Name.maxCardinality() << ')'<< std::endl; \
    if (TestInv( (Name), (seed))) { \
	std::cout << #Name << " failed !" << std::endl;	\
	return -1;					\
    }

    //TEST_INV(Modular<int8_t>, C17, 17); => GF(17) does not fit in Modular<int8_t>
    TEST_INV(Modular<int8_t>, C13, 13);
    TEST_INV(Modular<int16_t>, S17, 17);
    TEST_INV(Modular<int32_t>, Z17, 17);
    TEST_INV(Modular<int64_t>, LL17, 17);
    TEST_INV(Modular<uint8_t>, UC13, 13);
    TEST_INV(Modular<uint16_t>, US17, 17);
    TEST_INV(Modular<uint32_t>, UZ17, 17);
    TEST_INV(Modular<uint64_t>, ULL17, 17);
    TEST_INV(ModularCUS, CUS17, 17);
    TEST_INV(ModularSUZ, SUZ17, 17);
    TEST_INV(ModularZULL, ZULL17, 17);
    TEST_INV(ModularUCUS, UCUS17, 17);
    TEST_INV(ModularUSUZ, USUZ17, 17);
    TEST_INV(ModularUZULL, UZULL17, 17);
    TEST_INV(ModularFD, FD17, 17);
#ifdef __GIVARO_HAVE_INT128
    TEST_INV(ModularLLULLL, LLULLL17, 17);
    TEST_INV(ModularULLULLL, ULLULLL17, 17);
#endif

    TEST_INV(Modular<float>, F17, 17);
    TEST_INV(Modular<double>, D17, 17);
    TEST_INV(ModularExtended<float>, MEF17, 17);
    TEST_INV(ModularExtended<double>, MED17, 17);
    TEST_INV(Modular<Integer>, I17, 17);
    TEST_INV(Modular<RecInt::ruint128>, RU17, 17);
    TEST_INV(Modular<RecInt::rint128>, R17, 17);
    TEST_INV(Modular<Log16>, L17, 17);

    return 0;

}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
