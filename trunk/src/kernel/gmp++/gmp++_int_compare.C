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

/*! @file gmp++/gmp++_int_compare.C
 * @brief routines to compare integers.
 */

#ifndef __GIVARO_gmpxx_gmpxx_int_compare_C
#define __GIVARO_gmpxx_gmpxx_int_compare_C

#ifndef __GIVARO_INLINE_ALL
#include "gmp++/gmp++.h"
#endif
#include <cstdlib>

namespace Givaro {

	int compare(const Integer &a, const Integer& b)
	{
		return mpz_cmp ( (mpz_srcptr)&a.gmp_rep, (mpz_srcptr)&b.gmp_rep );
	}

	// absCompare
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

	int absCompare(const Integer &a, const uint64_t b)
	{
		return mpz_cmpabs_ui( (mpz_srcptr)&(a.gmp_rep), b);
	}

	int absCompare(const Integer &a, const uint32_t b)
	{
		return mpz_cmpabs_ui( (mpz_srcptr)&(a.gmp_rep), (uint64_t)b);
	}

	int absCompare(const Integer &a, const long int b)
	{
		return mpz_cmpabs_ui( (mpz_srcptr)&(a.gmp_rep), (uint64_t) std::abs(b));
	}

	int absCompare(const Integer &a, const int b)
	{
		return mpz_cmpabs_ui( (mpz_srcptr)&(a.gmp_rep), (uint64_t) std::abs(b));
	}

	// Operator !=
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

	int Integer::operator != (const uint32_t l) const
	{
		return mpz_cmp_ui((mpz_srcptr)&gmp_rep, (uint64_t) l) != 0;
	}

	int Integer::operator != (const long int l) const
	{
		return mpz_cmp_si((mpz_srcptr)&gmp_rep, l) != 0;
	}

	int Integer::operator != (const uint64_t l) const
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

	int operator != (long int l, const Integer& n)
	{
		return n.operator != (l);
	}

	int operator != (uint64_t l, const Integer& n)
	{
		return n.operator != (l);
	}

	int operator != (uint32_t l, const Integer& n)
	{
		return n.operator != (l);
	}

	// operator ==
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

	int Integer::operator == (const uint32_t l) const
	{
		return mpz_cmp_ui((mpz_srcptr)&gmp_rep, (uint64_t) l) == 0;
	}

	int Integer::operator == (const long int l) const
	{
		return mpz_cmp_si((mpz_srcptr)&gmp_rep, l) == 0;
	}

	int Integer::operator == (const uint64_t l) const
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

	int operator == (long int l, const Integer& n)
	{
		return n.operator == (l);
	}

	int operator == (uint64_t l, const Integer& n)
	{
		return n.operator == (l);
	}

	int operator == (uint32_t l, const Integer& n)
	{
		return n.operator == (l);
	}

	// Operator >
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

	int Integer::operator > (const uint32_t l) const
	{
		return mpz_cmp_ui((mpz_srcptr)&gmp_rep, (uint64_t) l) > 0;
	}

	int Integer::operator > (const long int l) const
	{
		return mpz_cmp_si((mpz_srcptr)&gmp_rep, l) > 0;
	}

	int Integer::operator > (const uint64_t l) const
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

	int operator > (long int l, const Integer& n)
	{
		return n.operator < (l);
	}

	int operator > (uint64_t l, const Integer& n)
	{
		return n.operator < (l);
	}

	int operator > (uint32_t l, const Integer& n)
	{
		return n.operator < (l);
	}

	// Operator <
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

	int Integer::operator < (const uint32_t l) const
	{
		return mpz_cmp_ui((mpz_srcptr)&gmp_rep, (uint64_t) l) < 0;
	}

	int Integer::operator < (const uint64_t l) const
	{
		return mpz_cmp_ui((mpz_srcptr)&gmp_rep, l) < 0;
	}

	int Integer::operator < (const int l) const
	{
		return mpz_cmp_si((mpz_srcptr)&gmp_rep, l) < 0;
	}

	int Integer::operator < (const long int l) const
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

	int operator < (long int l, const Integer& n)
	{
		return n.operator > (l);
	}

	int operator < (uint64_t l, const Integer& n)
	{
		return n.operator > (l);
	}

	int operator < (uint32_t l, const Integer& n)
	{
		return n.operator > (l);
	}

	// Operator >=
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

	int Integer::operator >= (const uint32_t l) const
	{
		return mpz_cmp_ui((mpz_srcptr)&gmp_rep, (uint64_t) l) >= 0;
	}

	int Integer::operator >= (const long int l) const
	{
		return mpz_cmp_si((mpz_srcptr)&gmp_rep, l) >= 0;
	}

	int Integer::operator >= (const uint64_t l) const
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

	int operator >= (long int l, const Integer& n)
	{
		return n.operator <= (l);
	}

	int operator >= (uint64_t l, const Integer& n)
	{
		return n.operator <= (l);
	}

	int operator >= (uint32_t l, const Integer& n)
	{
		return n.operator <= (l);
	}

	// Operator <=
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

	int Integer::operator <= (const uint32_t l) const
	{
		return mpz_cmp_ui((mpz_srcptr)&gmp_rep, (uint64_t) l) <= 0;
	}

	int Integer::operator <= (const uint64_t l) const
	{
		return mpz_cmp_ui((mpz_srcptr)&gmp_rep, l) <= 0;
	}

	int Integer::operator <= (const int l) const
	{
		return mpz_cmp_si((mpz_srcptr)&gmp_rep, l) <= 0;
	}

	int Integer::operator <= (const long int l) const
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

	int operator <= (long int l, const Integer& n)
	{
		return n.operator >= (l);
	}

	int operator <= (uint64_t l, const Integer& n)
	{
		return n.operator >= (l);
	}

	int operator <= (uint32_t l, const Integer& n)
	{
		return n.operator >= (l);
	}


	// compare to 1 and 0
	int isOne(const Integer& a)
	{
		return ! mpz_cmp_ui((mpz_srcptr)&(a.gmp_rep), 1UL);
	}
	int isMOne(const Integer& a)
	{
		return ! mpz_cmp_si((mpz_srcptr)&(a.gmp_rep), 1L);
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
	int isZero(const long int a)
	{
		return a ==0;
	}
	int isZero(const uint16_t a)
	{
		return a ==0;
	}
	int isZero(const uint32_t a)
	{
		return a ==0;
	}
	int isZero(const uint64_t a)
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
	int isZero(const long long int a)
	{
		return a ==0LL;
	}
#endif

}

#endif // __GIVARO_gmpxx_gmpxx_int_compare_C

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
