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

#ifndef __GIVARO_INLINE_ALL
#include "gmp++/gmp++.h"
#endif
#include <cstdlib>

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

	// absCompare
	/*! Compare the norm of two integers.
	 * @param a integer
	 * @param b integer
	 * @return \c 1 if \f$|a| > |b|\f$, \c 0 if \f$|a| = |b|\f$ and \p -1 otherwise.
	 */
	//@{
	int absCompare(const Integer &a, const Integer &b)
	{
		return mpz_cmpabs( (mpz_srcptr)&(a.gmp_rep), (mpz_srcptr)&(b.gmp_rep));
	}

	int absCompare(const Integer &a, const double b)
	{
		return mpz_cmpabs_d( (mpz_srcptr)&(a.gmp_rep), b);
	}

	int absCompare(const Integer &a, const float b)
	{
		return mpz_cmpabs_d( (mpz_srcptr)&(a.gmp_rep), (double)b);
	}

	int absCompare(const Integer &a, const unsigned long b)
	{
		return mpz_cmpabs_ui( (mpz_srcptr)&(a.gmp_rep), b);
	}

	int absCompare(const Integer &a, const unsigned b)
	{
		return mpz_cmpabs_ui( (mpz_srcptr)&(a.gmp_rep), (unsigned long)b);
	}

	int absCompare(const Integer &a, const long int b)
	{
		return mpz_cmpabs_ui( (mpz_srcptr)&(a.gmp_rep), (unsigned long) std::abs(b));
	}

	int absCompare(const Integer &a, const int b)
	{
		return mpz_cmpabs_ui( (mpz_srcptr)&(a.gmp_rep), (unsigned long) std::abs(b));
	}
	//@}

	// Operator !=

	/*! operator != (not equal)
	 * @param l integer
	 * @return \c 1 iff l == this
	 */
	//@{
	int Integer::operator != (const Integer & l) const
	{
		return mpz_cmp((mpz_srcptr)&gmp_rep,  (mpz_srcptr)l.get_mpz_const()) != 0;
	}

	int Integer::operator != (const double l) const
	{
		return mpz_cmp_d((mpz_srcptr)&gmp_rep,  l) != 0;
	}

	int Integer::operator != (const float l) const
	{
		return mpz_cmp_d((mpz_srcptr)&gmp_rep, (float)  l) != 0;
	}

	int Integer::operator != (const int l) const
	{
		return mpz_cmp_si((mpz_srcptr)&gmp_rep, l) != 0;
	}

	int Integer::operator != (const unsigned l) const
	{
		return mpz_cmp_ui((mpz_srcptr)&gmp_rep, (unsigned long) l) != 0;
	}

	int Integer::operator != (const long l) const
	{
		return mpz_cmp_si((mpz_srcptr)&gmp_rep, l) != 0;
	}

	int Integer::operator != (const unsigned long l) const
	{
		return mpz_cmp_ui((mpz_srcptr)&gmp_rep, l) != 0;
	}

	//-------------------------------------------------inline comparaison operators
	int operator != (double l, const Integer& n)
	{
		return n.operator != (l);
	}

	int operator != (float l, const Integer& n)
	{
		return n.operator != (l);
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

	int operator != (unsigned l, const Integer& n)
	{
		return n.operator != (l);
	}

	//@}

	//@{
	int Integer::operator == (const Integer & l) const
	{
		return mpz_cmp((mpz_srcptr)&gmp_rep,  (mpz_srcptr)l.get_mpz_const()) == 0;
	}

	int Integer::operator == (const double l) const
	{
		return mpz_cmp_d((mpz_srcptr)&gmp_rep,  l) == 0;
	}

	int Integer::operator == (const float l) const
	{
		return mpz_cmp_d((mpz_srcptr)&gmp_rep, (float)  l) == 0;
	}

	int Integer::operator == (const int l) const
	{
		return mpz_cmp_si((mpz_srcptr)&gmp_rep, l) == 0;
	}

	int Integer::operator == (const unsigned l) const
	{
		return mpz_cmp_ui((mpz_srcptr)&gmp_rep, (unsigned long) l) == 0;
	}

	int Integer::operator == (const long l) const
	{
		return mpz_cmp_si((mpz_srcptr)&gmp_rep, l) == 0;
	}

	int Integer::operator == (const unsigned long l) const
	{
		return mpz_cmp_ui((mpz_srcptr)&gmp_rep, l) == 0;
	}

	//-------------------------------------------------inline comparaison operators
	int operator == (double l, const Integer& n)
	{
		return n.operator == (l);
	}

	int operator == (float l, const Integer& n)
	{
		return n.operator == (l);
	}

	int operator == (int l, const Integer& n)
	{
		return n.operator == (l);
	}

	int operator == (long l, const Integer& n)
	{
		return n.operator == (l);
	}

	int operator == (unsigned long l, const Integer& n)
	{
		return n.operator == (l);
	}

	int operator == (unsigned l, const Integer& n)
	{
		return n.operator == (l);
	}

	//@}



	// Operator >
	//@{
	int Integer::operator > (const Integer & l) const
	{
		return mpz_cmp((mpz_srcptr)&gmp_rep,  (mpz_srcptr)l.get_mpz_const()) > 0;
	}

	int Integer::operator > (const double l) const
	{
		return mpz_cmp_d((mpz_srcptr)&gmp_rep,  l) > 0;
	}

	int Integer::operator > (const float l) const
	{
		return mpz_cmp_d((mpz_srcptr)&gmp_rep, (float)  l) > 0;
	}

	int Integer::operator > (const int l) const
	{
		return mpz_cmp_si((mpz_srcptr)&gmp_rep, l) > 0;
	}

	int Integer::operator > (const unsigned l) const
	{
		return mpz_cmp_ui((mpz_srcptr)&gmp_rep, (unsigned long) l) > 0;
	}

	int Integer::operator > (const long l) const
	{
		return mpz_cmp_si((mpz_srcptr)&gmp_rep, l) > 0;
	}

	int Integer::operator > (const unsigned long l) const
	{
		return mpz_cmp_ui((mpz_srcptr)&gmp_rep, l) > 0;
	}

	//-------------------------------------------------inline comparaison operators
	int operator > (double l, const Integer& n)
	{
		return n.operator < (l);
	}

	int operator > (float l, const Integer& n)
	{
		return n.operator < (l);
	}

	int operator > (int l, const Integer& n)
	{
		return n.operator < (l);
	}

	int operator > (long l, const Integer& n)
	{
		return n.operator < (l);
	}

	int operator > (unsigned long l, const Integer& n)
	{
		return n.operator < (l);
	}

	int operator > (unsigned l, const Integer& n)
	{
		return n.operator < (l);
	}

	// Operator <
	//@{
	int Integer::operator < (const Integer & l) const
	{
		return mpz_cmp((mpz_srcptr)&gmp_rep,  (mpz_srcptr)l.get_mpz_const()) < 0;
	}

	int Integer::operator < (const double l) const
	{
		return mpz_cmp_d((mpz_srcptr)&gmp_rep,  l) < 0;
	}

	int Integer::operator < (const float l) const
	{
		return mpz_cmp_d((mpz_srcptr)&gmp_rep, (float) l) < 0;
	}

	int Integer::operator < (const unsigned l) const
	{
		return mpz_cmp_ui((mpz_srcptr)&gmp_rep, (unsigned long) l) < 0;
	}

	int Integer::operator < (const unsigned long l) const
	{
		return mpz_cmp_ui((mpz_srcptr)&gmp_rep, l) < 0;
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
	int operator < (double l, const Integer& n)
	{
		return n.operator > (l);
	}

	int operator < (float l, const Integer& n)
	{
		return n.operator > (l);
	}

	int operator < (int l, const Integer& n)
	{
		return n.operator > (l);
	}

	int operator < (long l, const Integer& n)
	{
		return n.operator > (l);
	}

	int operator < (unsigned long l, const Integer& n)
	{
		return n.operator > (l);
	}

	int operator < (unsigned l, const Integer& n)
	{
		return n.operator > (l);
	}

	// Operator >=
	//@{
	int Integer::operator >= (const Integer & l) const
	{
		return mpz_cmp((mpz_srcptr)&gmp_rep,  (mpz_srcptr)l.get_mpz_const()) >= 0;
	}

	int Integer::operator >= (const double l) const
	{
		return mpz_cmp_d((mpz_srcptr)&gmp_rep,  l) >= 0;
	}

	int Integer::operator >= (const float l) const
	{
		return mpz_cmp_d((mpz_srcptr)&gmp_rep, (float)  l) >= 0;
	}

	int Integer::operator >= (const int l) const
	{
		return mpz_cmp_si((mpz_srcptr)&gmp_rep, l) >= 0;
	}

	int Integer::operator >= (const unsigned l) const
	{
		return mpz_cmp_ui((mpz_srcptr)&gmp_rep, (unsigned long) l) >= 0;
	}

	int Integer::operator >= (const long l) const
	{
		return mpz_cmp_si((mpz_srcptr)&gmp_rep, l) >= 0;
	}

	int Integer::operator >= (const unsigned long l) const
	{
		return mpz_cmp_ui((mpz_srcptr)&gmp_rep, l) >= 0;
	}

	//-------------------------------------------------inline comparaison operators
	int operator >= (double l, const Integer& n)
	{
		return n.operator <= (l);
	}

	int operator >= (float l, const Integer& n)
	{
		return n.operator <= (l);
	}

	int operator >= (int l, const Integer& n)
	{
		return n.operator <= (l);
	}

	int operator >= (long l, const Integer& n)
	{
		return n.operator <= (l);
	}

	int operator >= (unsigned long l, const Integer& n)
	{
		return n.operator <= (l);
	}

	int operator >= (unsigned l, const Integer& n)
	{
		return n.operator <= (l);
	}

	// Operator <=
	//@{
	int Integer::operator <= (const Integer & l) const
	{
		return mpz_cmp((mpz_srcptr)&gmp_rep,  (mpz_srcptr)l.get_mpz_const()) <= 0;
	}

	int Integer::operator <= (const double l) const
	{
		return mpz_cmp_d((mpz_srcptr)&gmp_rep,  l) <= 0;
	}

	int Integer::operator <= (const float l) const
	{
		return mpz_cmp_d((mpz_srcptr)&gmp_rep, (float) l) <= 0;
	}

	int Integer::operator <= (const unsigned l) const
	{
		return mpz_cmp_ui((mpz_srcptr)&gmp_rep, (unsigned long) l) <= 0;
	}

	int Integer::operator <= (const unsigned long l) const
	{
		return mpz_cmp_ui((mpz_srcptr)&gmp_rep, l) <= 0;
	}

	int Integer::operator <= (const int l) const
	{
		return mpz_cmp_si((mpz_srcptr)&gmp_rep, l) <= 0;
	}

	int Integer::operator <= (const long l) const
	{
		return mpz_cmp_si((mpz_srcptr)&gmp_rep, l) <= 0;
	}

	//-------------------------------------------------inline comparaison operators
	int operator <= (double l, const Integer& n)
	{
		return n.operator >= (l);
	}

	int operator <= (float l, const Integer& n)
	{
		return n.operator >= (l);
	}

	int operator <= (int l, const Integer& n)
	{
		return n.operator >= (l);
	}

	int operator <= (long l, const Integer& n)
	{
		return n.operator >= (l);
	}

	int operator <= (unsigned long l, const Integer& n)
	{
		return n.operator >= (l);
	}

	int operator <= (unsigned l, const Integer& n)
	{
		return n.operator >= (l);
	}



	// compare to 1 and 0
	int isOne(const Integer& a)
	{
		return ! mpz_cmp_ui((mpz_srcptr)&(a.gmp_rep), 1UL);
	}

	int nonZero(const Integer& a)
	{
		return mpz_cmp_ui((mpz_srcptr)&(a.gmp_rep), 0UL);
		// return (mpz_sgn((mpz_srcptr)&(a.gmp_rep)) != 0) ; // BB which one is faster ?
	}

	int isZero(const Integer& a)
	{
		return ! mpz_cmp_ui((mpz_srcptr)&(a.gmp_rep), 0UL);
		// return (mpz_sgn((mpz_srcptr)&(a.gmp_rep)) == 0) ; // BB which one is faster ?
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
