// Copyright(c)'2010 by The Givaro group
// This file is part of Givaro.
// written by BB
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

#include <iostream>

#include <givaro/givinteger.h>

//#define TEST_PASS
#undef TEST_PASS

#ifndef TEST_PASS
#define SONT_EQ(a,b)\
	if ( (a) != (b) ) { \
		std::cout << "erreur à la ligne " << __LINE__ << std::endl; \
		std::cout << a << "!=" << b << std::endl; \
		return -1 ; \
	}
#else
#define SONT_EQ(a,b)\
	if ( (a) != (b) ) { \
		std::cout << "erreur à la ligne " << __LINE__ << std::endl; \
		std::cout << a << "!=" << b << std::endl; \
	}
#endif

#define NB_ITERS 40

template<class U> inline bool IsNeg(const U p) { return (p<0); }
template<> inline bool IsNeg<unsigned long>(const unsigned long p) { return false; }

template< class T, class U>
int test1( const T m, const U p)
{/*{{{*/
	double pi (p);
	long int r = m % p;
	if (r<0)  r += (IsNeg(p) )?(-p):(p); // r est dans [[0,p-1]]
	const Integer M(m);
	const Integer P(p);
	Integer R ;

	Integer::mod(R,M,P);
	SONT_EQ(r,R);
	// R a bon r.

	Integer R1 = M%p ; //!XX
	SONT_EQ(R,R1);

	Integer R2 = M%P ;
	SONT_EQ(R,R2);

	Integer R3 = M%pi ; //!XX
	SONT_EQ(R,R3);

	R1 = M ;
	Integer::modin(R1,P);
	SONT_EQ(R,R1);


	R2 = M ;
	R2 %= P ;
	SONT_EQ(R,R2);

	return 0;
}/*}}}*/

int test2(Integer & M, Integer & P)
{/*{{{*/
	Integer RR ;
	Integer MM = M ;
	Integer PP = P ;
	mpz_mod(RR.get_mpz() ,MM.get_mpz() ,PP.get_mpz());
	if (RR<0) RR+=(PP<0)?(-PP):(PP);

	Integer R = 0 ;

	//!@todo existe pas !
	//R = Integer:: mod(M,P) ;
	Integer:: mod(R,M,P) ;
	SONT_EQ(RR,R);

	R = M ;
	Integer::modin(R,P);
	SONT_EQ(RR,R);

	R = M ;
	R %= P;
	SONT_EQ(RR,R);

	R = M ;
	R = M%P;
	SONT_EQ(RR,R);

	return 0;
}/*}}}*/


int main()
{/*{{{*/
#if (SIZEOF_LONG==8)
	long int m = 1253345363665346363;
#else
	long int m = 1665346363;
#endif

	long int p = 78678675;
	unsigned long int M(m);
	unsigned long int P(p);

	int rez = 0;

	rez =  test1(m,p);   if (rez) return 1 ;
	rez =  test1(-m,p);  if (rez) return 2 ;
	rez =  test1(m,-p);  if (rez) return 3 ;
	rez =  test1(-m,-p); if (rez) return 4 ;
	rez =  test1(m,P);   if (rez) return 5 ;
	rez =  test1(M,P);   if (rez) return 6 ;

	rez =  test1(m,m);   if (rez) return 1 ;
	rez =  test1(p,p);   if (rez) return 2 ;
	rez =  test1(-m,m);  if (rez) return 3 ;
	rez =  test1(-p,p);  if (rez) return 4 ;
	rez =  test1(m,-m);  if (rez) return 5 ;
	rez =  test1(p,-p);  if (rez) return 6 ;
	rez =  test1(-m,-m); if (rez) return 7 ;
	rez =  test1(-p,-p); if (rez) return 8 ;


	for (unsigned i = 0 ; i < NB_ITERS ; ++i)
	{/*{{{*/
		Integer M = Integer::random_between(680,700);
		Integer P = Integer::random_between(134,198);
		rez = test2(M,P); if (rez) return 7 ;
		Integer::negin(M);
		rez = test2(M,P); if (rez) return 8 ;
		Integer::negin(P);
		rez = test2(M,P); if (rez) return 9 ;
		Integer::negin(M);
		rez = test2(M,P); if (rez) return 10 ;
	}/*}}}*/
	//std::cout << "ok6" << std::endl;

	rez =  test1(0,p);   if (rez) return 11 ;
	rez =  test1(0,-p);   if (rez) return 12 ;
	//std::cout << "ok7" << std::endl;



	return 0;
}/*}}}*/

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen
