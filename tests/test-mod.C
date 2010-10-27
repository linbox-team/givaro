// Copyright(c)'2010 by The Givaro group
// This file is part of Givaro.
// written by BB
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software. 
// see the COPYRIGHT file for more details.

#include <iostream>

#include <givaro/givinteger.h>

#define SONT_EQ(a,b)\
	if ( (a) != (b) ) { \
		std::cout << "erreur à la ligne " << __LINE__ << std::endl; \
		std::cout << a << "!=" << b << std::endl; \
		return -1 ; \
	}

#define NB_ITERS 40

int main()
{/*{{{*/

	unsigned long int m = 1253345363665346363;
	unsigned long int p = 78678675;
	double pi (p);
	unsigned long int r = m % p;
	Integer M(m);
	Integer P(p);
	Integer R ;

	Integer::mod(R,M,P);
	SONT_EQ(R,r);

	Integer R1 = M%p ;
	SONT_EQ(R,R1);
	
	Integer R2 = M%P ;
	SONT_EQ(R,R2);

	Integer R3 = M%pi ;
	SONT_EQ(R,R3);

	for (unsigned i = 0 ; i < NB_ITERS ; ++i)
	{/*{{{*/
		M = Integer::random_between(680,700);
		P = Integer::random_between(134,198);

		//!@todo existe pas !
		//R = Integer:: mod(M,P) ; 
		Integer:: mod(R,M,P) ;

		Integer RR ;
		mpz_tdiv_r(RR.get_mpz() ,M.get_mpz() ,P.get_mpz());
		SONT_EQ(RR,R);

		Integer::negin(P);

		Integer:: mod(R,M,P) ;
		mpz_tdiv_r(RR.get_mpz() ,M.get_mpz() ,P.get_mpz());
		SONT_EQ(RR,R);

		Integer::negin(M);
		Integer:: mod(R,M,P) ;
		mpz_tdiv_r(RR.get_mpz() ,M.get_mpz() ,P.get_mpz());
		if (RR<0) RR-=P ;
		SONT_EQ(RR,R);

		Integer::negin(P);
		Integer:: mod(R,M,P) ;
		mpz_tdiv_r(RR.get_mpz() ,M.get_mpz() ,P.get_mpz());
		if (RR<0) RR+=P ;
		SONT_EQ(RR,R);
	}/*}}}*/


	return 0;
}/*}}}*/

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen
