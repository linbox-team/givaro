// Copyright(c)'1994-2025 by The Givaro group
// This file is part of Givaro.
// BB: adapted,enhanced from examples/FiniteFields/ff-arith.C
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

#include <iostream>

#include <givaro/modular.h>
#include <givaro/modular-balanced.h>
#include <givaro/montgomery.h>

#include <givaro/gfq.h>
#include <givaro/gfqext.h>
#include <givaro/extension.h>
#include <givaro/givintprime.h>

#include <recint/recint.h>

using namespace Givaro;

#define TESTE_EG( a, b )						\
    if (!F.areEqual((a),(b))) {						\
	F.write(F.write(std::cout,a) << "!=",b)				\
	    << " failed (at line " <<  __LINE__ << ")" << std::endl;	\
	return -1 ;							\
    }

#define JETESTE( a, seed )				\
    if (TestField( (a), int(seed)) ) {			\
	std::cout << #a << " failed !" << std::endl;	\
	return -1 ;					\
    }

#define JEONETESTE( F, a )				\
    if (TestOneField(F, (a))) {				\
	std::cout << #a << " failed !" << std::endl;	\
	return -1 ;					\
    }
    
template<class Field>
bool invertible(const Field& F, const typename Field::Element& a)
{
    auto ai(a);
    return F.mulin(F.inv(ai, a), a) == F.one;
}

template<class Field>
int TestOneField(const Field& F, const typename Field::Element& first)
{
#ifdef GIVARO_DEBUG
    F.write(std::cerr << "Testing ") << "(" << first << ") :" << std::endl;
#endif

    typename Field::Element a, b, c, d,a_,b_,c_,d_,ma;
    typename Field::Element e,e_;

    F.init(a, 0UL);
    TESTE_EG(a, F.zero);
    F.init(a, 1UL);
    //         F.write(std::cerr) << std::endl;
    //         F.write(std::cerr << "a: ", a) << std::endl;
    //         F.write(std::cerr << "1: ", F.one) << std::endl;
    TESTE_EG(a, F.one);

    F.inv(a_, a);

    TESTE_EG(a_, F.one);

    F.init(ma,-1L);
    //         F.write(std::cerr) << std::endl;
    //         F.write(std::cerr << "a: ", a) << std::endl;
    //         F.write(std::cerr << "ma: ", ma) << std::endl;
    //         F.write(std::cerr << "1: ", F.one) << std::endl;
    //         F.write(std::cerr << "-1: ", F.mOne) << std::endl;
    TESTE_EG(ma, F.mOne);

    F.inv(a_, ma);

    TESTE_EG(a_, F.mOne);

    F.assign(a, first);

    typename Field::RandIter g(F);
    while (!invertible(F, g.random(b)));

    F.init(c);            // empty constructor
    F.init(d);            // empty constructor

    F.add(c, a, b);       // c = a+b
    F.assign(c_,c);       // c_ <- c

    //         F.write(std::cerr) << std::endl;
    //         F.write(std::cerr << "a:=", a) << ';' << std::endl;
    //         F.write(std::cerr << "b:=", b) << ';' << std::endl;
    //         F.write(std::cerr << "c:=", c) << ';' << std::endl;
    //         F.write(std::cerr << "c_:=", c_) << ';' << std::endl;

    TESTE_EG(c,c_);
    F.subin(c_,a);
    //         F.write(std::cerr) << std::endl;
    //         F.write(std::cerr << "a:=", a) << ';' << std::endl;
    //         F.write(std::cerr << "b:=", b) << ';' << std::endl;
    //         F.write(std::cerr << "c:=", c) << ';' << std::endl;
    //         F.write(std::cerr << "c_:=", c_) << ';' << std::endl;

    TESTE_EG(b,c_);

    F.mul(c, a, b);     // c = a*b
    F.assign(c_,c);       // c_ <- c
    F.divin(c_,b);      // c_ == a ?

    //         F.write(std::cerr) << std::endl;
    //         F.write(std::cerr << "a: ", a) << std::endl;
    //         F.write(std::cerr << "b: ", b) << std::endl;
    //         F.write(std::cerr << "c: ", c) << std::endl;
    //         F.write(std::cerr << "c_: ", c_) << std::endl;
    TESTE_EG(a,c_);

    F.assign(c, a);
    F.mulin(c, b);     // c = a*b
    F.divin(c,b);      // c_ == a ?

    //         F.write(std::cerr) << std::endl;
    //         F.write(std::cerr << "a: ", a) << std::endl;
    //         F.write(std::cerr << "b: ", b) << std::endl;
    //         F.write(std::cerr << "c: ", c) << std::endl;
    TESTE_EG(a,c);

    F.axpy(d, a, b, c); // d = a*b + c;
    F.init(d_);
    F.axmy(d_,a,b,c); // d_ = a*b - c
    F.addin(d_,c);
    F.subin(d,c);

    //         F.write(std::cerr) << std::endl;
    //         F.write(std::cerr << "a:=", a) << ';' << std::endl;
    //         F.write(std::cerr << "b:=", b) << ';' << std::endl;
    //         F.write(std::cerr << "c:=", c) << ';' << std::endl;
    //         F.write(std::cerr << "d:=", d) << ';' << std::endl;
    //         F.write(std::cerr << "d_:=", d_) << ';' << std::endl;
    TESTE_EG(d_,d);

    F.sub(d,a,b); // d = a -b
    F.add(c,a,b); // c = a+b
    F.init(e);
    F.init(e_);
    F.mul(e,d,c); // e = d*c;
    F.mul(a_,a,a); // a_ = a*a
    F.mul(b_,b,b); // b_ = b*b
    F.sub(e_,a_,b_); // e_ = a_ - b_

    //         F.write(std::cerr) << std::endl;
    //         F.write(std::cerr << "a:=", a) << ';' << std::endl;
    //         F.write(std::cerr << "b:=", b) << ';' << std::endl;
    //         F.write(std::cerr << "c:=", c) << ';' << std::endl;
    //         F.write(std::cerr << "d:=", d) << ';' << std::endl;
    //         F.write(std::cerr << "e:=", e) << ';' << std::endl;
    //         F.write(std::cerr << "e_:=", e_) << ';' << std::endl;
    //         F.write(std::cerr << "a_:=", a_) << ';' << std::endl;
    //         F.write(std::cerr << "b_:=", b_) << ';' << std::endl;
    TESTE_EG(e,e_); // a^2 - b^2 = (a-b)(a+b) ;)

    // Four operations
    F.init(a_);
    F.assign(a_,a);
    F.addin(a, b) ;
    F.subin(a, b) ;
    F.mulin(a, b) ;
    F.divin(a, b) ;

    TESTE_EG(a_,a);

    F.maxpy(e, a, b, d); // e = d-a*b
    F.assign(e_,d);
    F.maxpyin(e_, a, b); // e = d - a*b;

    //         F.write(std::cerr << "a:=", a) << ';' << std::endl;
    //         F.write(std::cerr << "b:=", b) << ';' << std::endl;
    //         F.write(std::cerr << "d:=", d) << ';' << std::endl;
    //         F.write(std::cerr << "e:=", e) << ';' << std::endl;
    //         F.write(std::cerr << "e_:=", e_) << ';' << std::endl;
    TESTE_EG(e,e_);

    F.axmy(e, a, b, d); // e = a*b -d;
    F.assign(e_,d);
    F.maxpyin(e_, a, b); // e = d - a*b;
    F.negin(e_);

    TESTE_EG(e,e_);

#ifdef GIVARO_DEBUG
    F.write(std::cerr );
    std::cerr  << " done." << std::endl;
#endif
    return 0 ;

}

#ifndef NBITER
#define NBITER 50
#endif

template<class Field>
int TestField(const Field& F, const int seed)
{
    typename Field::Element x;
    typename Field::RandIter g(F, seed);
    
    F.init(x, 7UL);
    JEONETESTE(F,x);
    
    for (size_t i = 0; i< NBITER; ++i) {
	while (F.isZero(g.random(x)));
        JEONETESTE(F,x);
    }
    
    return 0;
}

template<class Ints>
Ints previousprime(const Ints& a) {
    static IntPrimeDom IPD;
    Integer aI(a);
    IPD.prevprimein(aI);
    return (Ints)(aI);
}

int main(int argc, char ** argv)
{
    int seed = int (argc>1?atoi(argv[1]):BaseTimer::seed());
#ifdef GIVARO_DEBUG
    std::cerr << "seed: " << seed << std::endl;
#endif
    Integer::seeding((unsigned long)seed);
    RecInt::srand(seed);
    
    using ModularCUS = Modular<int8_t, uint16_t>;
    using ModularSUZ = Modular<int16_t, uint32_t>;
    using ModularZULL = Modular<int32_t, uint64_t>;
    using ModularUCUS = Modular<uint8_t, uint16_t>;
    using ModularUSUZ = Modular<uint16_t, uint32_t>;
    using ModularUZULL = Modular<uint32_t, uint64_t>;
    // using ModularFD = Modular<float, double>;

#define TEST_SPECIFIC(Field, Name, Modulus...)		\
    Field Name(Modulus);				\
    JETESTE(Name, seed);

    //--------------------//
    //----- Modulo 2 -----//
    
    TEST_SPECIFIC(Modular<int8_t>, C2, 2);
    TEST_SPECIFIC(Modular<int16_t>, S2, 2);
    TEST_SPECIFIC(Modular<int32_t>, Z2, 2);
    TEST_SPECIFIC(Modular<int64_t>, LL2, 2);
    TEST_SPECIFIC(Modular<uint8_t>, UC2, 2);
    TEST_SPECIFIC(Modular<uint16_t>, US2, 2);
    TEST_SPECIFIC(Modular<uint32_t>, UZ2, 2);
    TEST_SPECIFIC(Modular<uint64_t>, ULL2, 2);
    TEST_SPECIFIC(ModularCUS, CUS2, 2);
    TEST_SPECIFIC(ModularSUZ, SUZ2, 2);
    TEST_SPECIFIC(ModularZULL, ZULL2, 2);
    TEST_SPECIFIC(ModularUCUS, UCUS2, 2);
    TEST_SPECIFIC(ModularUSUZ, USUZ2, 2);
    TEST_SPECIFIC(ModularUZULL, UZULL2, 2);
    TEST_SPECIFIC(Modular<Log16>, L2, 2);
    TEST_SPECIFIC(Modular<float>, F2, 2);
    TEST_SPECIFIC(Modular<double>, D2, 2);
    //TEST_SPECIFIC(ModularFD, FD2, 2);
    TEST_SPECIFIC(Modular<Integer>, I2, 2);
    TEST_SPECIFIC(Modular<RecInt::rint128>, R2, 2);
    TEST_SPECIFIC(Modular<RecInt::ruint128>, RU2, 2);

    //--------------------//
    //----- Modulo 3 -----//
    
    TEST_SPECIFIC(ModularBalanced<int32_t>, BZ3, 3);
    TEST_SPECIFIC(ModularBalanced<int64_t>, BLL3, 3);
    TEST_SPECIFIC(ModularBalanced<float>, BF3, 3);
    TEST_SPECIFIC(ModularBalanced<double>, BD3, 3);
    
    TEST_SPECIFIC(Montgomery<int32_t>, MZ3, 3);
    TEST_SPECIFIC(Montgomery<RecInt::ruint128>, MRU3, 3);

    //---------------------//
    //----- Modulo 13 -----//
    
    TEST_SPECIFIC(Modular<int8_t>, C13, 13);
    TEST_SPECIFIC(Modular<int16_t>, S13, 13);
    TEST_SPECIFIC(Modular<int32_t>, Z13, 13);
    TEST_SPECIFIC(Modular<int64_t>, LL13, 13);
    TEST_SPECIFIC(Modular<uint8_t>, UC13, 13);
    TEST_SPECIFIC(Modular<uint16_t>, US13, 13);
    TEST_SPECIFIC(Modular<uint32_t>, UZ13, 13);
    TEST_SPECIFIC(Modular<uint64_t>, ULL13, 13);
    TEST_SPECIFIC(ModularCUS, CUS13, 13);
    TEST_SPECIFIC(ModularSUZ, SUZ13, 13);
    TEST_SPECIFIC(ModularZULL, ZULL13, 13);
    TEST_SPECIFIC(ModularUCUS, UCUS13, 13);
    TEST_SPECIFIC(ModularUSUZ, USUZ13, 13);
    TEST_SPECIFIC(ModularUZULL, UZULL13, 13);
    TEST_SPECIFIC(Modular<Log16>, L13, 13);
    TEST_SPECIFIC(Modular<float>, F13, 13);
    TEST_SPECIFIC(Modular<double>, D13, 13);
    //TEST_SPECIFIC(ModularFD, FD13, 13);
    TEST_SPECIFIC(Modular<Integer>, I13, 13);
    TEST_SPECIFIC(Modular<RecInt::rint128>, R13, 13);
    TEST_SPECIFIC(Modular<RecInt::ruint128>, RU13, 13);
    
    TEST_SPECIFIC(ModularBalanced<int32_t>, BZ13, 13);
    TEST_SPECIFIC(ModularBalanced<int64_t>, BLL13, 13);
    TEST_SPECIFIC(ModularBalanced<float>, BF13, 13);
    TEST_SPECIFIC(ModularBalanced<double>, BD13, 13);
    
    TEST_SPECIFIC(Montgomery<int32_t>, MZ13, 13);
    TEST_SPECIFIC(Montgomery<RecInt::ruint128>, MRU13, 13);

    //--------------------------------//
    //----- Modulo maximal prime -----//

#define TEST_LAST_PRIME(Field, Name)			\
    Field Name(previousprime(Field::getMaxModulus()));	\
    JETESTE(Name, seed);
    
    TEST_LAST_PRIME(Modular<int8_t>, Cpmax);
    TEST_LAST_PRIME(Modular<int16_t>, Spmax);
    TEST_LAST_PRIME(Modular<int32_t>, Zpmax);
    TEST_LAST_PRIME(Modular<int64_t>, LLpmax);
    TEST_LAST_PRIME(Modular<uint8_t>, UCpmax);
    TEST_LAST_PRIME(Modular<uint16_t>, USpmax);
    TEST_LAST_PRIME(Modular<uint32_t>, UZpmax);
    TEST_LAST_PRIME(Modular<uint64_t>, ULLpmax);
    TEST_LAST_PRIME(ModularCUS, CUSpmax);
    TEST_LAST_PRIME(ModularSUZ, SUZpmax);
    TEST_LAST_PRIME(ModularZULL, ZULLpmax);
    TEST_LAST_PRIME(ModularUCUS, UCUSpmax);
    TEST_LAST_PRIME(ModularUSUZ, USUZpmax);
    TEST_LAST_PRIME(ModularUZULL, UZULLpmax);
    TEST_LAST_PRIME(Modular<Log16>, Lpmax);
    TEST_LAST_PRIME(Modular<float>, Fpmax);
    TEST_LAST_PRIME(Modular<double>, Dpmax);
    //TEST_LAST_PRIME(ModularFD, FDpmax);
    TEST_LAST_PRIME(Modular<RecInt::rint128>, Rpmax);
    TEST_LAST_PRIME(Modular<RecInt::ruint128>, RUpmax);
    
    TEST_LAST_PRIME(ModularBalanced<int32_t>, BZpmax);
    TEST_LAST_PRIME(ModularBalanced<int64_t>, BLLpmax);
    TEST_LAST_PRIME(ModularBalanced<float>, BFpmax);
    TEST_LAST_PRIME(ModularBalanced<double>, BDpmax);
    
    TEST_LAST_PRIME(Montgomery<int32_t>, MZpmax);
    TEST_LAST_PRIME(Montgomery<RecInt::ruint128>, MRUpmax);

    //-------------------------//
    //----- Galois fields -----//
    
    TEST_SPECIFIC(GFqDom<int32_t>, GF13, 13);
    TEST_SPECIFIC(GFqDom<int32_t>, GFpmax, 65521UL);
    TEST_SPECIFIC(GFqDom<int64_t>, GFLLpmax, 4194301ULL);

    // Zech log finite field with 256 elements
    // and prescribed 1 + x +x^3 +x^4 +x^8 irreducible polynomial
    std::vector< GFqDom<int64_t>::Residu_t > Irred(9);
    Irred[0] = 1; Irred[1] = 1; Irred[2] = 0; Irred[3] = 1;
    Irred[4] = 1; Irred[5] = 0; Irred[6] = 0; Irred[7] = 0;
    Irred[8] = 1;
    TEST_SPECIFIC(GFqDom<int64_t>, GF256, 2, 8, Irred);

    TEST_SPECIFIC(GFqDom<int32_t>, GF625, 5, 4);
    TEST_SPECIFIC(GFqExt<int32_t>, GF81, 3, 4);

    // Zech log finite field with 2Mb tables
    TEST_SPECIFIC(GFqDom<int64_t>, GF2M, 2, 20);
    TEST_SPECIFIC(GFqDom<int64_t>, GF2M1, 2, 2);
    TEST_SPECIFIC(GFqDom<int64_t>, GF11E3, 11, 3);
    TEST_SPECIFIC(Extension<GFqDom<int64_t>>, GF11E9, GF11E3, 3);
    TEST_SPECIFIC(Extension<>, GF13E8, 13, 8);

    return 0;
}

