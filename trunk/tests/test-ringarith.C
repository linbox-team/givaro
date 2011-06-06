// Copyright(c)'1994-2009,2010 by The Givaro group
// This file is part of Givaro.
// BB: adapted,enhanced from examples/FiniteFields/ff-arith.C
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

#include <iostream>
#include "givaro/givzpz16std.h"
#include "givaro/givzpz16table1.h"
#include "givaro/givzpz32std.h"
#include "givaro/givzpz32uns.h"
#include <givaro/givzpzInt.h>
#include <givaro/givzpz64std.h>
#include <givaro/givzpz.h>
#include <givaro/givpoly1.h>

#include <givaro/givinteger.h>
using namespace Givaro;

#ifdef GIVARO_DEBUG
long long TTcount = 0;
#endif


#define TESTE_EG( a, b ) \
if (!F.areEqual((a),(b))) {\
	F.write( F.write(std::cout,a) << "!=",b) << " failed (at line " <<  __LINE__ << ")" << std::endl; \
	return(-1); \
}

#define JETESTE( a, s ) \
if (TestRing( (a), (s)) ) {\
	std::cout << #a << " failed !" << std::endl;\
	return -1 ; \
}

#define JEPOLTESTE( a, s ) \
if (TestPolRing( (a), (s) ) ) {\
	std::cout << #a << " failed !" << std::endl;\
	return -1 ; \
}

#define JEONETESTE( F, a, x ) \
if (TestOneRing(F,a,x)) {\
	std::cout << #a << " failed !" << std::endl;\
	return -1 ; \
}


template<class R, class T1, class T2>
struct InitOrAssign {
    void operator()(const R& r, T1& t1, const T2& t2) {
        r.init(t1,t2);
    }
};


template<class R, class T>
struct InitOrAssign<R,T,T> {
    void operator()(const R& r, T& t1, const T& t2) {
        r.assign(t1,t2);
    }
};




template<class Ring, class T1, class T2>
int TestOneRing(const Ring& F, const T1 FIRSTINT, const T2 FIRSTFLOAT)
{/*{{{*/
#ifdef GIVARO_DEBUG
	std::cerr << "testing " ;
	F.write(std::cerr );
	std::cerr  << " : " << std::flush;
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

	InitOrAssign<Ring,typename Ring::Element,T1>()(F, a, FIRSTINT);
        InitOrAssign<Ring,typename Ring::Element,T2>()(F, b, FIRSTFLOAT);
//       F.init(a, FIRSTINT);
// 	F.init(b, FIRSTFLOAT);

	F.init(c);            // empty constructor
	F.init(d);            // empty constructor

	F.add(c, a, b);       // c = a+b
	F.init(c_);           //! @warning F.init(c_,c); ne marche pas !
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

	TESTE_EG(e,e_) // a^2 - b^2 = (a-b)(a+b) ;



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
	F.axmyin(e_, a, b); // e = a*b-e=a*b-d;

	F.negin(e_);

	TESTE_EG(e,e_);



#ifdef GIVARO_DEBUG
	F.write(std::cerr );
	std::cerr  << " done." << std::endl;
        ++TTcount;
#endif
	return 0 ;

}/*}}}*/

#define NBITER 50

template<class Ring>
int TestRing(const Ring& F, const int seed)
{/*{{{*/
    long ch = (long) F.characteristic();
    JEONETESTE(F,7UL,-29.3);
    srand48(seed);
    for(size_t i=0; i< NBITER; ++i) {
        typename Ring::Element x;
        float d;
       	do {
		d = float((double)ch*drand48()) ;
		F.init(x,d );
        } while(F.isZero(x));
        int a; do {
            F.init(x, a = (int)lrand48());
        } while(F.isZero(x));
        JEONETESTE(F,x,d);
    }
    return 0;
}/*}}}*/

#define DEGMAX 75
#define NBITERD 10

template<class Ring>
int TestPolRing(const Ring& F, const int seed)
{/*{{{*/
    GivRandom generator(seed);
    srand48(seed);

    for(size_t i=0; i< NBITERD; ++i) {
        int d1 = int (lrand48() % DEGMAX);
        int d2 = int (lrand48() % DEGMAX);
#ifdef GIVARO_DEBUG
        std::cout << d1 << ' ' << d2 << ' ';
#endif
        typename Ring::Element x, d;
        do {
            F.random(generator, x, Degree(d1));
        } while(F.isZero(x));
        do {
            F.random(generator, d, Degree(d2));
        } while(F.isZero(x));
        JEONETESTE(F,x,d);
    }
    return 0;
}/*}}}*/

int main(int argc, char ** argv)
{/*{{{*/
    int seed = int(argc>1?atoi(argv[1]):BaseTimer::seed());
#ifdef GIVARO_DEBUG
    std::cerr << "seed: " << seed << std::endl;
#endif
    Integer::seeding(seed);


	// modulo 13 over 16 bits
	ZpzDom<Std16> C13(13);
	JETESTE(C13,seed);

	// modulo 13 over 32 bits
	ZpzDom<Std32> Z13(13);
	JETESTE(Z13,seed);

	// modulo 13 over unsigned 32 bits
	ZpzDom<Unsigned32> U13(13);
	JETESTE(U13,seed);

#ifdef __USE_Givaro_SIXTYFOUR__
	// modulo 13 over 64 bits
	ZpzDom<Std64> LL13(13UL);
	JETESTE(LL13,seed);
#endif

	// modulo 13 over arbitrary size
	ZpzDom<Integer> IntZ13(13);
	JETESTE(IntZ13,seed);


	// modulo 2 over 16 bits
	ZpzDom<Std16> C2(2);
	JETESTE(C2,seed);

	// modulo 2 over 32 bits
	ZpzDom<Std32> Z2(2);
	JETESTE(Z2,seed);

	// modulo 2 over unsigned 32 bits
	ZpzDom<Unsigned32> U2(2);
	JETESTE(U2,seed);

#ifdef __USE_Givaro_SIXTYFOUR__
	// modulo 2 over 64 bits
	ZpzDom<Std64> LL2(2UL);
	JETESTE(LL2,seed);
#endif

	// modulo 2 fully tabulated
	ZpzDom<Log16> L2(2);
	JETESTE(L2,seed);

	// modulo 13 over arbitrary size
	ZpzDom<Integer> IntZ2(2);
	JETESTE(IntZ2,seed);

// --------------------------------------------
	// modulo 4 over 16 bits
	ZpzDom<Std16> C4(4);
	JETESTE(C4,seed);

	// modulo 4 over 32 bits
	ZpzDom<Std32> Z4(4);
	JETESTE(Z4,seed);

	// modulo 4 over unsigned 32 bits
	ZpzDom<Unsigned32> U4(4);
	JETESTE(U4,seed);

#ifdef __USE_Givaro_SIXTYFOUR__
	// modulo 2 over 64 bits
	ZpzDom<Std64> LL4(4UL);
	JETESTE(LL4,seed);
#endif

	// modulo 13 over arbitrary size
	ZpzDom<Integer> IntZ4(4);
	JETESTE(IntZ4,seed);

// --------------------------------------------
	// modulo 75 over 16 bits
	ZpzDom<Std16> C75(75);
	JETESTE(C75,seed);

	// modulo 75 over 32 bits
	ZpzDom<Std32> Z75(75);
	JETESTE(Z75,seed);

	// modulo 75 over unsigned 32 bits
	ZpzDom<Unsigned32> U75(75);
	JETESTE(U75,seed);

#ifdef __USE_Givaro_SIXTYFOUR__
	// modulo 2 over 675 bits
	ZpzDom<Std64> LL75(75UL);
	JETESTE(LL75,seed);
#endif

	// modulo 13 over arbitrary size
	ZpzDom<Integer> IntZ75(75);
	JETESTE(IntZ75,seed);


        Poly1Dom< ZpzDom<Std16>, Dense > CP13(C13, "X");
	JETESTE(CP13,seed); JEPOLTESTE(CP13,seed);
        Poly1Dom< ZpzDom<Std32>, Dense > ZP13(Z13, "X");
	JETESTE(ZP13,seed); JEPOLTESTE(ZP13,seed);

        Poly1Dom< ZpzDom<Unsigned32>, Dense > UP13(U13, "X");
	JETESTE(UP13,seed); JEPOLTESTE(UP13,seed);

        Poly1Dom< ZpzDom<Std64>, Dense > LLP13(LL13, "X");
	JETESTE(LLP13,seed); JEPOLTESTE(LLP13,seed);

        Poly1Dom< ZpzDom<Integer>, Dense > IntZP13(IntZ13, "X");
	JETESTE(IntZP13,seed); JEPOLTESTE(IntZP13,seed);



        Poly1Dom< ZpzDom<Std16>, Dense > CP75(C75, "X");
	JEPOLTESTE(CP75,seed);

        Poly1Dom< ZpzDom<Std32>, Dense > ZP75(Z75, "X");
	JEPOLTESTE(ZP75,seed);

        Poly1Dom< ZpzDom<Unsigned32>, Dense > UP75(U75, "X");
	JEPOLTESTE(UP75,seed);

        Poly1Dom< ZpzDom<Std64>, Dense > LLP75(LL75, "X");
	JEPOLTESTE(LLP75,seed);

        Poly1Dom< ZpzDom<Integer>, Dense > IntZP75(IntZ75, "X");
	JEPOLTESTE(IntZP75,seed);



        Poly1Dom< Poly1Dom< ZpzDom<Integer>, Dense >, Dense> IntZPP75(IntZP75, "Y");
	JEPOLTESTE(IntZPP75,seed);

        Poly1Dom< Poly1Dom< Poly1Dom< ZpzDom<Integer>, Dense >, Dense>, Dense > IntZPPP75(IntZPP75, "Z");
	JEPOLTESTE(IntZPPP75,seed);

#ifdef GIVARO_DEBUG
        std::cerr << std::endl << "Success:" << TTcount << std::endl;
#endif

	return 0;
}/*}}}*/

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen
