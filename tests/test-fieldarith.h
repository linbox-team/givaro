// Copyright(c)'1994-2025 by The Givaro group
// This file is part of Givaro.
// BB: adapted,enhanced from examples/FiniteFields/ff-arith.C
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

#include <iostream>

#include <givaro/givinteger.h>
#include <givaro/givintprime.h>

using namespace Givaro;

#define TEST_EQUALITY( a, b )						\
if (!F.areEqual((a),(b))) {						\
    F.write(F.write(std::cerr,a) << "!=",b)				\
    << " failed (at line " <<  __LINE__ << ")" << std::endl;	\
    return -1 ;							\
}

#define TEST_FIELD_SEVERAL_TIMES_SEEDED( a, seed )				\
if (TestField( (a), int(seed)) ) {			\
    std::cerr << #a << " failed !" << std::endl;	\
    return -1 ;					\
}

#define TEST_ONE_FIELD_SEEDED( F, a )				\
if (TestOneField(F, (a))) {				\
    std::cerr << #a << " failed !" << std::endl;	\
    return -1 ;					\
} else std::clog << "PASSED." << std::endl;

#define TEST_LAST_PRIME(Field, Name)			\
Field Name(previousprime(Field::maxCardinality()));	\
TEST_FIELD_SEVERAL_TIMES_SEEDED(Name, seed);

#define TEST_SPECIFIC(Field, Name, Modulus...)		\
Field Name(Modulus);				\
TEST_FIELD_SEVERAL_TIMES_SEEDED(Name, seed);


template<class Field>
bool invertible(const Field& F, const typename Field::Element& a)
{
    //     auto ai(a);
    //     return F.mulin(F.inv(ai, a), a) == F.one;
    Integer ai;
    F.convert(ai,a);
    return (gcd(ai,Integer(F.characteristic()))==1);
}

template<class Field>
int TestOneField(const Field& F, const typename Field::Element& first)
{
    F.write(std::clog << "Testing ");
    F.write(std::clog<< "(", first) << ") :";

    typename Field::Element a, b, c, d,a_,b_,c_,d_,ma;
    typename Field::Element e,e_;

    F.init(a, 0);
    TEST_EQUALITY(a, F.zero);

    F.init(a, 0u);
    TEST_EQUALITY(a, F.zero);

    F.init(a, 1);
    //         F.write(std::cerr) << std::endl;
    //         F.write(std::cerr << "a: ", a) << std::endl;
    //         F.write(std::cerr << "1: ", F.one) << std::endl;
    TEST_EQUALITY(a, F.one);

    F.init(a, 1u);
    TEST_EQUALITY(a, F.one);

    F.inv(a_, a);
    TEST_EQUALITY(a_, F.one);

    F.init(ma,1); F.negin(ma);
    //         F.write(std::cerr) << std::endl;
    //         F.write(std::cerr << "a: ", a) << std::endl;
    //         F.write(std::cerr << "ma: ", ma) << std::endl;
    //         F.write(std::cerr << "1: ", F.one) << std::endl;
    //         F.write(std::cerr << "-1: ", F.mOne) << std::endl;
    TEST_EQUALITY(ma, F.mOne);

    F.init(ma,1_i64); F.negin(ma);
    TEST_EQUALITY(ma, F.mOne);

    F.inv(a_, ma);

    TEST_EQUALITY(a_, F.mOne);

#ifdef __GIVARO_DEBUG
    // F.write(std::cerr) << std::endl;
    // F.write(std::cerr << "0: ", F.zero) << std::endl;

    try {
        F.inv(a, F.zero);
    } catch(const GivMathDivZero& e) {
        //        std::cerr << "Correctly catched division by zero: " << e << std::endl;
    }
    catch (...) {
        F.mulin(a, F.zero);
        TEST_EQUALITY(a, F.one);
        std::cerr << "Error division by zero allowed even in DEBUG MODE" << std::endl;
    }


    // F.write(std::cerr << "1/0: ", a) << std::endl;
#endif


    F.assign(a, first);


    typename Field::RandIter g(F);
    while (!invertible(F, g.random(b))) {}

    F.init(c);            // empty constructor
    F.init(d);            // empty constructor

    F.add(c, a, b);       // c = a+b
    F.assign(c_,c);       // c_ <- c


    TEST_EQUALITY(c,c_);
    F.subin(c_,a);
    //         F.write(std::cerr) << std::endl;
    //         F.write(std::cerr << "a:=", a) << ';' << std::endl;
    //         F.write(std::cerr << "b:=", b) << ';' << std::endl;
    //         F.write(std::cerr << "c:=", c) << ';' << std::endl;
    //         F.write(std::cerr << "c_:=", c_) << ';' << std::endl;

    TEST_EQUALITY(b,c_);

    F.mul(c, a, b);     // c = a*b
    F.assign(c_,c);       // c_ <- c
    F.divin(c_,b);      // c_ == a ?

    //        F.write(std::cerr) << std::endl;
    //        F.write(std::cerr << "a: ", a) << std::endl;
    //        F.write(std::cerr << "b: ", b) << std::endl;
    //        F.write(std::cerr << "c: ", c) << std::endl;
    //        F.write(std::cerr << "c_: ", c_) << std::endl;
    TEST_EQUALITY(a,c_);

    F.assign(c, a);
    F.mulin(c, b);     // c = a*b
    F.divin(c,b);      // c_ == a ?

    //         F.write(std::cerr) << std::endl;
    //         F.write(std::cerr << "a: ", a) << std::endl;
    //         F.write(std::cerr << "b: ", b) << std::endl;
    //         F.write(std::cerr << "c: ", c) << std::endl;
    TEST_EQUALITY(a,c);

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
    TEST_EQUALITY(d_,d);

    F.sub(d,a,b); // d = a -b
    F.add(c,a,b); // c = a+b
    F.init(e);
    F.init(e_);
    F.mul(e,d,c); // e = d*c;
    F.mul(a_,a,a); // a_ = a*a
    F.mul(b_,b,b); // b_ = b*b
    F.sub(e_,a_,b_); // e_ = a_ - b_
    // F.write(std::cerr) << std::endl;
    // F.write(std::cerr << "a:=", a) << ';' << std::endl;
    // F.write(std::cerr << "b:=", b) << ';' << std::endl;
    // F.write(std::cerr << "c:=", c) << ';' << std::endl;
    // F.write(std::cerr << "d:=", d) << ';' << std::endl;
    // F.write(std::cerr << "e:=", e) << ';' << std::endl;
    // F.write(std::cerr << "e_:=", e_) << ';' << std::endl;
    // F.write(std::cerr << "a_:=", a_) << ';' << std::endl;
    // F.write(std::cerr << "b_:=", b_) << ';' << std::endl;
    TEST_EQUALITY(e,e_); // a^2 - b^2 = (a-b)(a+b) ;)

    // Four operations
    F.init(a_);
    F.assign(a_,a);
    F.addin(a, b) ;
    F.subin(a, b) ;
    F.mulin(a, b) ;
    F.divin(a, b) ;

    TEST_EQUALITY(a_,a);

    F.maxpy(e, a, b, d); // e = d-a*b
    F.assign(e_,d);
    F.maxpyin(e_, a, b); // e = d - a*b;

    //         F.write(std::cerr << "a:=", a) << ';' << std::endl;
    //         F.write(std::cerr << "b:=", b) << ';' << std::endl;
    //         F.write(std::cerr << "d:=", d) << ';' << std::endl;
    //         F.write(std::cerr << "e:=", e) << ';' << std::endl;
    //         F.write(std::cerr << "e_:=", e_) << ';' << std::endl;
    TEST_EQUALITY(e,e_);

    F.axmy(e, a, b, d); // e = a*b -d;
    F.assign(e_,d);
    F.maxpyin(e_, a, b); // e = d - a*b;
    F.negin(e_);

    TEST_EQUALITY(e,e_);

#ifdef __GIVARO_DEBUG
    F.write(std::cerr );
    std::cerr  << " done." << std::endl;
#endif
    return 0 ;

}

#ifndef NBITER
#define NBITER 50
#endif

template<class Field>
int TestField(const Field& F, const uint64_t seed)
{
    typename Field::Element x;
    typename Field::RandIter g(F, seed);

    F.init(x, 1);
    TEST_ONE_FIELD_SEEDED(F,x);

    for (size_t i = 0; i< NBITER; ++i) {
        while (F.isZero(g.random(x))) {}
        TEST_ONE_FIELD_SEEDED(F,x);
    }

    return 0;
}

template<class Ints>
Ints previousprime(const Ints& a) {
    static Givaro::IntPrimeDom IPD;
    Integer aI(a);
    IPD.prevprimein(aI);
    return (Ints)(aI);
}



/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
