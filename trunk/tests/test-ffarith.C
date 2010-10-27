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

#define TESTE_EG( a, b ) \
if (!F.areEqual((a),(b))) {\
	std::cout << F.write(std::cout,a) << "!=" << F.write(std::cout,b) << " failed (at line " <<  __LINE__ << ")" << std::endl; \
	return(-1); \
}

#define JETESTE( a ) \
if (TestField( (a) )) {\
	std::cout << #a << " failed !" << std::endl;\
	return -1 ; \
}

template<class Field>
int TestField(const Field& F) 
{/*{{{*/
	std::cerr << "testing " ; 
	F.write(std::cerr );
	std::cerr  << " : " << std::flush;

	typename Field::Element a, b, c, d,a_,b_,c_,d_;
	typename Field::Element e,e_;
	F.init(a, 7UL);
	F.init(b, -29.3);

	F.init(c);            // empty constructor
	F.init(d);            // empty constructor

	F.add(c, a, b);       // c = a+b
	F.init(c_);           //! @warning F.init(c_,c); ne marche pas !
	F.assign(c_,c);       // c_ <- c
	TESTE_EG(c,c_);
	F.subin(c_,a);
	TESTE_EG(b,c_);

	F.mul(c, a, b);     // c = a*b
	F.assign(c_,c);       // c_ <- c
	F.divin(c_,b);      // c_ == b ?
	TESTE_EG(a,c_);

	F.axpy(d, a, b, c); // d = a*b + c;
	F.init(d_);
	F.axmy(d_,a,b,c); // d_ = a*b - c
	F.addin(d_,c);
	F.subin(d,c);

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

	TESTE_EG(e,e_);


	F.axmy(e, a, b, d); // e = a*b -d;

	F.assign(e_,d);
	F.axmyin(e_, a, b); // e = d - a*b;

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

	return 0 ;

}/*}}}*/


int main(int argc, char ** argv) 
{/*{{{*/

	// modulo 13 over 16 bits
	ZpzDom<Std16> C13(13); 
	JETESTE(C13);

	// modulo 13 over 32 bits
	ZpzDom<Std32> Z13(13); 
	JETESTE(Z13);

	// modulo 13 over unsigned 32 bits
	ZpzDom<Unsigned32> U13(13); 
	JETESTE(U13);

#ifdef __USE_Givaro_SIXTYFOUR__
	// modulo 13 over 64 bits
	ZpzDom<Std64> LL13(13UL); 
	JETESTE(LL13);
#endif

	// modulo 13 fully tabulated
	ZpzDom<Log16> L13(13); 
	JETESTE(L13);

	// modulo 13 over 32 bits with Montgomery reduction
	Montgomery<Std32> M13(13); 
	JETESTE(M13);

	// modulo 2 over 16 bits
	ZpzDom<Std16> C2(2); 
	JETESTE(C2);

	// modulo 2 over 32 bits
	ZpzDom<Std32> Z2(2);
	JETESTE(Z2);

	// modulo 2 over unsigned 32 bits
	ZpzDom<Unsigned32> U2(2); 
	JETESTE(U2);

#ifdef __USE_Givaro_SIXTYFOUR__
	// modulo 2 over 64 bits
	ZpzDom<Std64> LL2(2UL); 
	JETESTE(LL2);
#endif

	// modulo 2 fully tabulated
	ZpzDom<Log16> L2(2); 
	JETESTE(L2);

	// modulo 2 over 32 bits with Montgomery reduction
	Montgomery<Std32> M2(2); 
	JETESTE(M2);

	Montgomery<Std32> M3(39989);
	JETESTE(M3);

	// modulo 13 with primitive root representation
	GFqDom<int> GF13( 13 ); 
	JETESTE(GF13);

	// modulo 13 over arbitrary size
	ZpzDom<Integer> IntZ13(13);
	JETESTE(IntZ13);

	// Zech log finite field with 5^4 elements
	GFqDom<int> GF625( 5, 4 );
	JETESTE(GF625);

	// Zech log finite field with 3^4 elements
	// Using the Q-adic Transform
	GFqExt<int> GF81( 3, 4 );
	JETESTE(GF81);

	// Zech log finite field with 2Mb tables
#ifndef __GIVARO__DONOTUSE_longlong__
	GFqDom<long long> GF2M( 2, 20 );
#else
	GFqDom<long> GF2M( 2, 20) ;
#endif
	JETESTE(GF2M);

#ifndef __GIVARO__DONOTUSE_longlong__
	GFqDom<long long> GF2M1( 2, 2 );
#else
	GFqDom<long> GF2M1( 2, 2) ;
#endif
	JETESTE(GF2M1);


	std::cerr << std::endl ;


	return 0;
}/*}}}*/

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen
