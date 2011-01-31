// Copyright(c)'1994-2009,2010 by The Givaro group
// This file is part of Givaro.
// BB: adapted,enhanced from examples/FiniteFields/ff-arith.C
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

#include <iostream>
#include <givaro/givzpz64std.h>
#include <givaro/givzpz.h>
#include <givaro/givgfq.h>
#include <givaro/givmontg32.h>
#include <givaro/givgfqext.h>
#include <givaro/givextension.h>

#define TESTE_EG( a, b ) \
if (!F.areEqual((a),(b))) {\
	std::cout << F.write(std::cout,a) << "!=" << F.write(std::cout,b) << " failed (at line " <<  __LINE__ << ")" << std::endl; \
	return(-1); \
}

#define JETESTE( a, s ) \
if (TestField( (a), (s)) ) {\
	std::cout << #a << " failed !" << std::endl;\
	return -1 ; \
}

#define JEONETESTE( F, a, x ) \
if (TestOneField(F,a,x)) {\
	std::cout << #a << " failed !" << std::endl;\
	return -1 ; \
}

template<class Field>
int TestOneField(const Field& F, const int FIRSTINT, const float FIRSTFLOAT)
{/*{{{*/
#ifdef GIVARO_DEBUG
	std::cerr << "testing " ;
	F.write(std::cerr );
	std::cerr  << " : " << std::flush;
#endif



	typename Field::Element a, b, c, d,a_,b_,c_,d_;
	typename Field::Element e,e_;

        F.init(a, 0UL);
        TESTE_EG(a, F.zero);
        F.init(a, 1UL);
//         F.write(std::cerr) << std::endl;
//         F.write(std::cerr << "a: ", a) << std::endl;
//         F.write(std::cerr << "1: ", F.one) << std::endl;
        TESTE_EG(a, F.one);

	F.init(a, FIRSTINT);
	F.init(b, FIRSTFLOAT);

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

	F.mul(c, a, b);     // c = a*b
	F.assign(c_,c);       // c_ <- c
	F.divin(c_,b);      // c_ == a ?

//         F.write(std::cerr) << std::endl;
//         F.write(std::cerr << "a: ", a) << std::endl;
//         F.write(std::cerr << "b: ", b) << std::endl;
//         F.write(std::cerr << "c: ", c) << std::endl;
//         F.write(std::cerr << "c_: ", c_) << std::endl;
	TESTE_EG(a,c_);

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

	TESTE_EG(e,e_) // a^2 - b^2 = (a-b)(a+b) ;)

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


	if ((unsigned) F.characteristic() != (unsigned)2) {
		F.init(e,1);
		F.init(a,22996);
		F.inv(b,a);
		F.mul(c,b,a);

		TESTE_EG(e,c);

		F.init(a,22996);
		F.init(b,22996);
		F.invin(a);
		F.mulin(a,b);
	}

	TESTE_EG(e,a);

	F.init(a,37403);
	F.inv(b,a);
	F.mul(c,b,a);

	TESTE_EG(e,c);

	F.init(a,37403);
	F.init(b,37403);
	F.invin(a);
	F.mulin(a,b);

	TESTE_EG(e,a);

#ifdef GIVARO_DEBUG
	F.write(std::cerr );
	std::cerr  << " done." << std::endl;
#endif
	return 0 ;

}/*}}}*/

#define NBITER 50

template<class Field>
int TestField(const Field& F, const int seed)
{/*{{{*/
    long ch = F.characteristic();
    JEONETESTE(F,7UL,-29.3);
    srand48(seed);
    for(size_t i=0; i< NBITER; ++i) {
        typename Field::Element x;
        float d; do {
            F.init(x, d = ch*drand48());
        } while(F.isZero(x));
        int a; do {
            F.init(x, a = lrand48());
        } while(F.isZero(x));
        JEONETESTE(F,a,d);
    }
    return 0;
}/*}}}*/



int main(int argc, char ** argv)
{/*{{{*/
    int seed = (argc>1?atoi(argv[1]):BaseTimer::seed());
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

	// modulo 13 fully tabulated
	ZpzDom<Log16> L13(13);
	JETESTE(L13,seed);

	// modulo 13 over 32 bits with Montgomery reduction
	Montgomery<Std32> M13(13);
	JETESTE(M13,seed);


// Maximal values

	// prime modulo max over 16 bits
	ZpzDom<Std16> Cmax(251);
	JETESTE(Cmax,seed);

	// prime modulo max fully tabulated
	ZpzDom<Log16> Lmax(251);
	JETESTE(Lmax,seed);

	// prime modulo max over 32 bits
	ZpzDom<Std32> Zmax(65521);
	JETESTE(Zmax,seed);

	// prime modulo max over 32 bits
	ZpzDom<Unsigned32> Umax(65521);
	JETESTE(Umax,seed);

	// prime modulo max over 32 bits with Montgomery reduction
	Montgomery<Std32> Mmax(40499);
	JETESTE(Mmax,seed);

#ifdef __USE_Givaro_SIXTYFOUR__
	// prime modulo max over 64 bits
	ZpzDom<Std64> LLmax(4294967291ULL);
	JETESTE(LLmax,seed);
#endif

      


// Characteristic 2


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

	// modulo 3 over 32 bits with Montgomery reduction
	Montgomery<Std32> M2(3);
	JETESTE(M2,seed);

	Montgomery<Std32> M3(39989);
	JETESTE(M3,seed);

	// modulo 13 with primitive root representation
	GFqDom<int> GF13( 13 );
	JETESTE(GF13,seed);

	// modulo 13 over arbitrary size
	ZpzDom<Integer> IntZ13(13);
	JETESTE(IntZ13,seed);

	// Zech log finite field with 5^4 elements
	GFqDom<int> GF625( 5, 4 );
	JETESTE(GF625,seed);

	// Zech log finite field with 3^4 elements
	// Using the Q-adic Transform
	GFqExt<int> GF81( 3, 4 );
	JETESTE(GF81,seed);

	// Zech log finite field with 2Mb tables
#ifndef __GIVARO__DONOTUSE_longlong__
	GFqDom<long long> GF2M( 2, 20 );
#else
	GFqDom<long> GF2M( 2, 20) ;
#endif
	JETESTE(GF2M,seed);

#ifndef __GIVARO__DONOTUSE_longlong__
	GFqDom<long long> GF2M1( 2, 2 );
#else
	GFqDom<long> GF2M1( 2, 2) ;
#endif
	JETESTE(GF2M1,seed);


        Extension<> GF13E8(13,8);
	JETESTE(GF13E8,seed);

#ifndef __GIVARO__DONOTUSE_longlong__
        GFqDom<long long> GF11e3( 11, 3 );
	JETESTE(GF11e3,seed);
        Extension<GFqDom<long long> > GF11e9(GF11e3,3);
	JETESTE(GF11e9,seed);
#else
        GFqDom<long> GF11e3( 11, 3) ;
	JETESTE(GF11e3,seed);
        Extension<> GF11e9(GF11e3,3);
	JETESTE(GF11e9,seed);
#endif

#ifdef GIVARO_DEBUG
	std::cerr << std::endl ;
#endif


	return 0;
}/*}}}*/

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen
