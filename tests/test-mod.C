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

template< class T, class U>
int test3( const T m, const U p)
{/*{{{*/
	int pi (p);
	long int q = m / p;
	const Integer M(m);
	const Integer P(p);
	Integer Q ;

	Integer::div(Q,M,P);
	SONT_EQ(q,Q);

	Integer Q1 = M/p ; //!XX
	SONT_EQ(Q,Q1);

	Integer Q2 = M/P ;
	SONT_EQ(Q,Q2);

	Integer Q3 = M/pi ; //!XX
	SONT_EQ(Q,Q3);

	Q1 = M ;
	Integer::divin(Q1,P);
	SONT_EQ(Q,Q1);


	Q2 = M ;
	Q2 /= P ;
	SONT_EQ(Q,Q2);

	return 0;
}/*}}}*/

#include <cassert>

int main()
{/*{{{*/
#if (SIZEOF_LONG==8)
	long int m = 1253345363665346363;
#else
	long int m = 1665346363;
#endif

	long int p = 78678675;
	std::cout << "mod " << p << std::endl;
	unsigned long int M(m);
	unsigned long int P(p);

        Integer mone(-1);
            // CONDITION: mpz_tdiv_ui does NOT consider the sign of gmp_rep
        assert(mpz_tdiv_ui( (mpz_ptr)&mone, 3) == 1);

	for (long i = -7 ; i < 7 ; ++i){
		long j = 3 ;
		std::cout << i/j << '=';
		Integer I = i;
		Integer::divin(I,j);
		std::cout << I << '|' ;

		std::cout << i/(-j) << '=';
		I = i;
		Integer::divin(I,-j);
		std::cout << I << std::endl ;


	}


	int rez = 0;

	rez =  test1(m,p);   if (rez) return 1 ;
	rez =  test1(-m,p);  if (rez) return 2 ;
	rez =  test1(m,-p);  if (rez) return 3 ;
	rez =  test1(-m,-p); if (rez) return 4 ;

	rez =  test1(m,P);   if (rez) return 5 ;
	rez =  test1(m,-P);  if (rez) return 6 ;
	rez =  test1(-m,P);  if (rez) return 7 ;
	rez =  test1(-m,-P); if (rez) return 8 ;

	rez =  test1(M,P);   if (rez) return 9 ;
	rez =  test1(M,-P);  if (rez) return 10 ;
	rez =  test1(-M,P);  if (rez) return 11 ;
	rez =  test1(-M,-P); if (rez) return 12 ;

	rez =  test1(-M,p);  if (rez) return 13 ;
	rez =  test1(M,-p);  if (rez) return 14 ;
	rez =  test1(M,p);   if (rez) return 15 ;
	rez =  test1(-M,-p); if (rez) return 16 ;

	rez =  test1(m,m);   if (rez) return 16 ;
	rez =  test1(-m,m);  if (rez) return 18 ;
	rez =  test1(m,-m);  if (rez) return 19 ;
	rez =  test1(-m,-m); if (rez) return 20 ;

	rez =  test1(P,P);   if (rez) return 21 ;
	rez =  test1(P,-P);  if (rez) return 22 ;
	rez =  test1(-P,P);  if (rez) return 23 ;
	rez =  test1(-P,-P); if (rez) return 24 ;


	rez =  test3(m,p);   if (rez) return 11 ;
	rez =  test3(-m,p);  if (rez) return 12 ;
	rez =  test3(m,-p);  if (rez) return 13 ;
	rez =  test3(-m,-p); if (rez) return 14 ;

	rez =  test3(m,P);   if (rez) return 15 ;
	rez =  test3(m,-P);  if (rez) return 16 ;
	rez =  test3(-m,P);  if (rez) return 17 ;
	rez =  test3(-m,-P); if (rez) return 18 ;

	rez =  test3(M,P);   if (rez) return 19 ;
	rez =  test3(M,-P);  if (rez) return 110 ;
	rez =  test3(-M,P);  if (rez) return 111 ;
	rez =  test3(-M,-P); if (rez) return 112 ;

	rez =  test3(-M,p);  if (rez) return 113 ;
	rez =  test3(M,-p);  if (rez) return 114 ;
	rez =  test3(M,p);   if (rez) return 115 ;
	rez =  test3(-M,-p); if (rez) return 116 ;

	rez =  test3(m,m);   if (rez) return 116 ;
	rez =  test3(-m,m);  if (rez) return 118 ;
	rez =  test3(m,-m);  if (rez) return 119 ;
	rez =  test3(-m,-m); if (rez) return 120 ;

	rez =  test3(P,P);   if (rez) return 121 ;
	rez =  test3(P,-P);  if (rez) return 122 ;
	rez =  test3(-P,P);  if (rez) return 123 ;
	rez =  test3(-P,-P); if (rez) return 124 ;



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

	rez =  test1(0,p);    if (rez) return 21 ;
	rez =  test1(0,-p);   if (rez) return 22 ;

	rez =  test3(0,p);    if (rez) return 23 ;
	rez =  test3(0,-p);   if (rez) return 24 ;

	//std::cout << "ok7" << std::endl;



	return 0;
}/*}}}*/

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen
