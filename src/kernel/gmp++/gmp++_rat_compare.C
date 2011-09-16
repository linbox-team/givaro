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

#ifndef __GIVARO_gmpxx_gmpxx_rat_compare_C
#define __GIVARO_gmpxx_gmpxx_rat_compare_C

#include "gmp++/gmp++.h"

#if 0
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

}
#endif

namespace Givaro {

		Rationel abs(const Rationel& f)
		{
			Rationel absf ;
			mpq_abs((mpq_ptr)absf.get_mpq(),(mpq_srcptr)f.get_mpq());
			return absf ;
		}
}

#endif // __GIVARO_gmpxx_gmpxx_rat_compare_C

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
