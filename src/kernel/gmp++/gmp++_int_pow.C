// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int_pow.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: JG Dumas
// $Id: gmp++_int_pow.C,v 1.4 2009-09-17 14:28:22 jgdumas Exp $
// ==========================================================================
// Description:
/** @file gmp++/gmp++_int_pow.C
 * powing stuff.
 */

#ifndef __GIVARO_gmpxx_gmpxx_int_pow_C
#define __GIVARO_gmpxx_gmpxx_int_pow_C

#ifndef __GIVARO_INLINE_ALL
#include "gmp++/gmp++.h"
#endif
#include <cstdlib>

namespace Givaro {

	int isperfectpower(const Integer& n)
	{
		return mpz_perfect_power_p((mpz_srcptr)&(n.gmp_rep));
	}

	Integer& pow(Integer& Res, const uint64_t n, const uint64_t p)
	{
		mpz_ui_pow_ui( (mpz_ptr)&(Res.gmp_rep), n, p);
		return Res;
	}
	Integer& pow(Integer& Res, const Integer& n, const uint64_t p)
	{
		__gmpz_pow_ui( (mpz_ptr)&(Res.gmp_rep), (mpz_srcptr)&n.gmp_rep, p);
		return Res;
	}

	Integer pow(const Integer& n, const uint64_t p)
	{
		if (p == 0) return Integer::one;

		Integer Res;
		return pow(Res,n,p);
	}

	Integer& pow(Integer& Res, const Integer& n, const long int l)
	{
		// Beware of negative values
		return pow(Res, n, (uint64_t) std::abs(l) );
	}
	Integer pow(const Integer& n, const long int l)
	{
		if (l < 0)
			return Integer::zero;
		return pow(n, (uint64_t) std::abs(l) );
	}

	Integer& powmod(Integer& Res, const Integer& n, const uint64_t p, const Integer& m)
	{
		mpz_powm_ui( (mpz_ptr)&(Res.gmp_rep), (mpz_srcptr)&n.gmp_rep, p, (mpz_srcptr)&m.gmp_rep);
		return Res;
	}

	Integer powmod(const Integer& n, const uint64_t p, const Integer& m)
	{
		if (p == 0) return Integer::one;
		Integer Res;
		return powmod(Res,n,p,m);
	}

	Integer& powmod(Integer& Res, const Integer& n, const long e, const Integer& m)
	{
		if (e < 0) {
			inv(Res, n, m);
			return powmod(Res, Res, (uint64_t)std::abs(e), m);
		}
		else {
			return powmod (Res, n, (uint64_t)(e), m);
		}
	}
	Integer powmod(const Integer& n, const long e, const Integer& m)
	{
		Integer Res;
		return powmod(Res, n, e, m);
	}


	Integer& powmod(Integer& Res, const Integer& n, const Integer& e, const Integer& m)
	{
		mpz_powm( (mpz_ptr)&(Res.gmp_rep), (mpz_srcptr)&n.gmp_rep, (mpz_srcptr)&e.gmp_rep, (mpz_srcptr)&m.gmp_rep);
		return Res;
	}
	Integer powmod(const Integer& n, const Integer& e, const Integer& m)
	{
		if (e == 0) return Integer::one;
		if (e < 0)  return Integer::zero;
		Integer Res;
		return powmod(Res, n, e, m);
	}

}
#endif // __GIVARO_gmpxx_gmpxx_int_pow_C
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen
