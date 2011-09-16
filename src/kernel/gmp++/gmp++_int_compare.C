// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int_compare.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: M. Samama, T. Gautier
// $Id: gmp++_int_compare.C,v 1.6 2009-09-17 14:28:22 jgdumas Exp $
// ==========================================================================

/*! @file givaro/gmp++_int_compare.C
 * @brief routines to compare integers.
 */

#ifndef __GIVARO_gmpxx_gmpxx_int_compare_C
#define __GIVARO_gmpxx_gmpxx_int_compare_C

#include "gmp++/gmp++.h"

namespace Givaro {

	/*! Compares two integers.
	 * @param a integer
	 * @param b integer
	 * @return \c 1 if \f$a > b\f$, \c 0 if \f$a = b\f$ and \p -1 otherwise.
	 */
	int compare(const Integer &a, const Integer& b)
	{
		return mpz_cmp ( (mpz_srcptr)&a.gmp_rep, (mpz_srcptr)&b.gmp_rep );
	}

	/*! Compare the norm of two integers.
	 * @param a integer
	 * @param b integer
	 * @return \c 1 if \f$|a| > |b|\f$, \c 0 if \f$|a| = |b|\f$ and \p -1 otherwise.
	 */
	int absCompare(const Integer &a, const Integer &b)
	{
		return mpz_cmpabs( (mpz_srcptr)&(a.gmp_rep), (mpz_srcptr)&(b.gmp_rep));
	}

	int Integer::operator != (const int l) const
	{
		return mpz_cmp_si ( (mpz_srcptr)&gmp_rep, l ) != 0;
	}

	int Integer::operator != (const long l) const
	{
		return mpz_cmp_si ( (mpz_srcptr)&gmp_rep, l ) != 0;
	}

	//unsigned long ops added by Dan Roche, 6-26-04
	int Integer::operator != (const unsigned long l) const
	{
		return mpz_cmp_ui ( (mpz_srcptr)&gmp_rep, l ) != 0;
	}

	int Integer::operator > (const unsigned long l) const
	{
		return mpz_cmp_ui((mpz_srcptr)&gmp_rep, l) > 0;
	}

	int Integer::operator < (const unsigned long l) const
	{
		return mpz_cmp_ui((mpz_srcptr)&gmp_rep, l) < 0;
	}

	int Integer::operator > (const int l) const
	{
		return mpz_cmp_si((mpz_srcptr)&gmp_rep, l) > 0;
	}

	int Integer::operator > (const long l) const
	{
		return mpz_cmp_si((mpz_srcptr)&gmp_rep, l) > 0;
	}

	int Integer::operator < (const int l) const
	{
		return mpz_cmp_si((mpz_srcptr)&gmp_rep, l) < 0;
	}

	int Integer::operator < (const long l) const
	{
		return mpz_cmp_si((mpz_srcptr)&gmp_rep, l) < 0;
	}

	//-------------------------------------------------inline comparaison operators
	int operator != (const Integer& a , const Integer& b)
	{
		return compare(a,b) != 0;
	}

	 int operator != (int l, const Integer& n)
	{
		return n.operator != (l);
	}

	 int operator != (long l, const Integer& n)
	{
		return n.operator != (l);
	}

	 int operator != (unsigned long l, const Integer& n)
	{
		return n.operator != (l);
	}

	 int operator == (const Integer& a, const Integer& b)
	{
		return compare(a,b) == 0;
	}

	 int operator == (int l, const Integer& n)
	{
		return (! (n.operator != (l)));
	}

	 int operator == (long l, const Integer& n)
	{
		return (! (n.operator != (l)));
	}

	 int operator == (unsigned long l, const Integer& n)
	{
		return (! (n.operator != (l)));
	}

	 int operator == (const Integer& n, unsigned long l)
	{
		return (! (n.operator != (l)));
	}

	 int operator == (const Integer& n, int l)
	{
		return (! (n.operator != (l)));
	}

	 int operator == (const Integer& n, long l)
	{
		return (! (n.operator != (l)));
	}

	 int operator < (const Integer& a , const Integer& b)
	{
		return compare(a,b) < 0;
	}

	 int operator < (const int l, const Integer& n)
	{
		return n > l;
	}

	 int operator < (const long l, const Integer& n)
	{
		return n > l;
	}

	 int operator < (const unsigned long l, const Integer& n)
	{
		return n > l;
	}

	 int operator <= (const Integer& n, unsigned long l)
	{
		return (! (n > l) );
	}

	 int operator <= (unsigned long l, const Integer& n)
	{
		return (! (n < l) );
	}

	 int operator >= (unsigned long l, const Integer& n)
	{
		return (! (n < l) );
	}

	 int operator >= (const Integer& n, unsigned long l)
	{
		return (! (n < l) );
	}

	 int operator > (int l, const Integer& n)
	{
		return n < l;
	}

	 int operator > (long l, const Integer& n)
	{
		return n < l;
	}

	 int operator > (unsigned long l, const Integer& n)
	{
		return n < l;
	}

	 int operator >  (const Integer& a , const Integer& b)
	{
		return compare(a,b) > 0;
	}

	 int operator <= (const Integer& a, const Integer& b)
	{
		return compare(a,b) <= 0;
	}

	 int operator <= (const Integer& n, int l)
	{
		return (! (n > l) );
	}

	 int operator <= (const Integer& n, long l)
	{
		return (! (n > l) );
	}

	 int operator <= (int l, const Integer& n)
	{
		return (! (n < l) );
	}

	 int operator <= (long l, const Integer& n)
	{
		return (! (n < l) );
	}

	 int operator >= (const Integer& a, const Integer& b)
	{
		return compare(a,b) >= 0;
	}

	 int operator >= (int l, const Integer& n)
	{
		return (! (n > l) );
	}

	 int operator >= (long l, const Integer& n)
	{
		return (! (n > l) );
	}

	 int operator >= (const Integer& n, int l)
	{
		return (! (n < l) );
	}

	 int operator >= (const Integer& n, long l)
	{
		return (! (n < l) );
	}

	 // compare to 1 and 0
	  int isOne(const Integer& a)
	{
		return ! mpz_cmp_ui((mpz_srcptr)&(a.gmp_rep), 1UL);
	}

	 int isZero(const Integer& a)
	{
		return ! mpz_cmp_ui((mpz_srcptr)&(a.gmp_rep), 0UL);
		// return (mpz_sgn((mpz_srcptr)&(a.gmp_rep)) == 0) ; // BB which one is faster ?
	}

	 int nonZero(const Integer& a)
	{
		return mpz_cmp_ui((mpz_srcptr)&(a.gmp_rep), 0UL);
		// return (mpz_sgn((mpz_srcptr)&(a.gmp_rep)) != 0) ; // BB which one is faster ?
	}

	 int isZero(const short int a)
	{
		return a ==0;
	}
	 int isZero(const int a)
	{
		return a ==0;
	}
	 int isZero(const long a)
	{
		return a ==0;
	}
	 int isZero(const unsigned short int a)
	{
		return a ==0;
	}
	 int isZero(const unsigned int a)
	{
		return a ==0;
	}
	 int isZero(const unsigned long a)
	{
		return a ==0UL;
	}
#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
#if 1 /*  use of C++0x long long integer constant */
	 int isZero(const unsigned long long a)
	{
		return a ==0ULL;
	}
#endif
	 int isZero(const long long a)
	{
		return a ==0LL;
	}
#endif

}

#endif // __GIVARO_gmpxx_gmpxx_int_compare_C

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
