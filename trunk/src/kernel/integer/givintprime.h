// =================================================================== //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <06 Jun 06 14:48:16 Jean-Guillaume.Dumas@imag.fr>
// =================================================================== //


/*! @file givintprime.h
 * @ingroup integers
 * @brief primes
 * - Prime numbers
 * - Modular powering,
 * - Fermat numbers,
 * - Primality tests
 * - Factorization : (There are parameters to fix)
 * .
 */
#ifndef __GIVARO_integers_prime_H
#define __GIVARO_integers_prime_H

#ifndef _GIVARO_ISPRIMETESTS_
#define _GIVARO_ISPRIMETESTS_ 5
#endif

#include "givaro/givinteger.h"

namespace Givaro {

	// =================================================================== //
	//! Fermat numbers
	// =================================================================== //
	class FermatDom : public IntegerDom {
	public:
		FermatDom() : IntegerDom() {}
		Rep& fermat (Rep&, const long)  const ;
		int pepin (const long) const ;
	};


	// =================================================================== //
	// Primality tests and factorization algorithms
	// =================================================================== //

	// Those macros are parameters to fix

	// primes known
	// first array
#define LOGMAX 3512
#define TABMAX 32768
	// second array
#define LOGMAX2 3031
#define TABMAX2 65536
	// Bounds between big and small
#define BOUNDARY_isprime TABMAX
#define BOUNDARY_2_isprime TABMAX2

	// =================================================================== //
	//! Primality tests
	// =================================================================== //
	class IntPrimeDom : public IntegerDom {
	public:
		IntPrimeDom() :  IntegerDom() {}

		int isprime(const Rep& n, int r=_GIVARO_ISPRIMETESTS_) const
		{
			/*
			   return probab_prime(n);
			   */
			//             return ((n)<BOUNDARY_isprime ?  isprime_Tabule(n) :
			//                     (n)<BOUNDARY_2_isprime ? isprime_Tabule2(n) :
			//                     probab_prime(n));
			long l;
			return int (int(islt(n,BOUNDARY_isprime) ?  isprime_Tabule((int)convert(l,n)):
					islt(n,BOUNDARY_2_isprime) ? isprime_Tabule2((int)convert(l,n)):
					local_prime(n,r)));
		}

		// if p is a prime power, p = r^return
		// else return is 0 and r is undefined
		unsigned int isprimepower(Rep&, const Rep&) const ;

		template<class RandIter>
		unsigned int Miller(RandIter& g, const Rep& n=_GIVARO_ISPRIMETESTS_) const  ;

		template<class RandIter>
		Rep& test_Lehmann(RandIter& g, Rep&, const Rep& n=_GIVARO_ISPRIMETESTS_) const  ;

		template<class RandIter>
		int Lehmann(RandIter& g, const Rep& n=_GIVARO_ISPRIMETESTS_)  const ;

		int isprime_Tabule(const int n) const ;
		int isprime_Tabule2(const int n) const ;

		Rep& nextprime(Rep&, const Rep&, int r=_GIVARO_ISPRIMETESTS_) const ;
		Rep& prevprime(Rep&, const Rep&, int r=_GIVARO_ISPRIMETESTS_) const ;
		Rep& nextprimein(Rep&, int r=_GIVARO_ISPRIMETESTS_) const ;
		Rep& prevprimein(Rep&, int r=_GIVARO_ISPRIMETESTS_) const ;


		// Using Integer
		int local_prime(const Rep& n, int r=_GIVARO_ISPRIMETESTS_) const
		{
			return probab_prime(n,r);
		}

	private:
		static int IP[LOGMAX+5];  // -- table for Tabule
		static const int * TP;    // -- shifted table
		static int IP2[LOGMAX2+5]; // -- table for Tabule2
		static const int * TP2;    // -- shifted table
#if 0
		   static int Tabule2(const Integer& p) ;
		   static int Tabule(const Integer& p) ;
		   static int _memTab2[LOGMAX2+5];   // -- table for Tabule2
		   static const int* _Tab2; // -- shifted _memTabule2
		   static int _memTab[];    // -- table for Tabule
		   static const int* _Tab;  // -- shifted _memTabule
#endif
	};

} // Givaro
#include "givaro/givintprime.inl"
#endif // __GIVARO_integers_prime_H

/*  -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
