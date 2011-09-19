// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int_sub.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: M. Samama, T. Gautier
// $Id: gmp++_int_sub.C,v 1.4 2009-09-17 14:28:22 jgdumas Exp $
// ==========================================================================


#ifndef __GIVARO_gmpxx_gmpxx_int_sub_C
#define __GIVARO_gmpxx_gmpxx_int_sub_C

#include "gmp++/gmp++.h"

namespace Givaro {

	//-------------------------------------------------- operator -
	Integer& Integer::subin(Integer& res, const Integer& n)
	{
		if (isZero(n)) return res;
		if (isZero(res)) return res = - n;
		mpz_sub( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&n.gmp_rep );
		return res;
	}
	Integer& Integer::subin(Integer& res, const long n)
	{
		if (isZero(n)) return res;
		if (isZero(res)) return res = - n;
		int sgn = Givaro::sign(n);
		if (sgn >0) mpz_sub_ui( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, n);
		else mpz_add_ui((mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, -n);
		return res;
	}
	Integer& Integer::subin(Integer& res, const unsigned long n)
	{
		if (isZero(n)) return res;
		if (isZero(res)) return res = - n;
		mpz_sub_ui( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&res.gmp_rep, n);
		return res;
	}

	Integer& Integer::sub(Integer& res, const Integer& n1, const Integer& n2)
	{
		if (isZero(n1)) return res = - n2;
		if (isZero(n2)) return res = n1;
		mpz_sub( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&n1.gmp_rep, (mpz_srcptr)&n2.gmp_rep);
		return res;
	}
	Integer& Integer::sub(Integer& res, const Integer& n1, const long n2)
	{
		if (isZero(n1)) return res = - n2;
		if (isZero(n2)) return res = n1;
		int sgn = Givaro::sign(n2);
		if (sgn >0) mpz_sub_ui( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&n1.gmp_rep, n2);
		else mpz_add_ui((mpz_ptr)&res.gmp_rep, (mpz_srcptr)&n1.gmp_rep, -n2);
		return res;
	}
	Integer& Integer::sub(Integer& res, const Integer& n1, const unsigned long n2)
	{
		if (isZero(n1)) return res = - n2;
		if (isZero(n2)) return res = n1;
		mpz_sub_ui( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&n1.gmp_rep, n2);
		return res;
	}

	Integer& Integer::neg(Integer& res, const Integer& n)
	{
		mpz_neg( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&n.gmp_rep);
		return res;
	}

	Integer& Integer::negin(Integer& res)
	{
		mpz_neg( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep);
		return res;
	}


	Integer& Integer::operator -= (const Integer& n)
	{
		if (isZero(n)) return *this;
		if (isZero(*this)) return logcpy(-n);
		//   Rep (res.gmp_rep)( MAX(SZ_REP(n.gmp_rep),SZ_REP(gmp_rep)) );
		mpz_sub( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, (mpz_srcptr)&n.gmp_rep) ;
		return *this;
	}

	Integer& Integer::operator -= (const unsigned long l)
	{
		if (l==0) return *this;
		if (isZero(*this)) return logcpy(Integer(-l));
		//   Rep (res.gmp_rep)( MAX(SZ_REP(gmp_rep),1) );
		mpz_sub_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, l);
		return *this;
	}

	Integer& Integer::operator -= (const long l)
	{
		if (l==0) return *this;
		if (isZero(*this)) return logcpy(Integer(-l));
		//   Rep (res.gmp_rep)( MAX(SZ_REP(gmp_rep),1) );
		int sgn = Givaro::sign(l);
		if (sgn >0) mpz_sub_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, l);
		else mpz_add_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, -l);
		return *this;
	}


	Integer Integer::operator - (const Integer& n) const
	{
		if (isZero(n)) return *this;
		if (isZero(*this)) return -n;
		//   Rep (res.gmp_rep)( MAX(SZ_REP(n.gmp_rep),SZ_REP(gmp_rep)) );
		Integer res;
		mpz_sub( (mpz_ptr)&(res.gmp_rep), (mpz_srcptr)&gmp_rep, (mpz_srcptr)&n.gmp_rep) ;
		return res;
	}

	Integer Integer::operator - (const unsigned long l) const
	{
		if (l==0) return *this;
		if (isZero(*this)) return Integer(-l);
		//   Rep (res.gmp_rep)( MAX(SZ_REP(gmp_rep),1) );
		Integer res;
		mpz_sub_ui( (mpz_ptr)&(res.gmp_rep), (mpz_srcptr)&gmp_rep, l);
		return res;
	}

	Integer Integer::operator - (const long l) const
	{
		if (l==0) return *this;
		if (isZero(*this)) return Integer(-l);
		//   Rep (res.gmp_rep)( MAX(SZ_REP(gmp_rep),1) );
		Integer res;
		int sgn = Givaro::sign(l);
		if (sgn >0) mpz_sub_ui( (mpz_ptr)&(res.gmp_rep), (mpz_srcptr)&gmp_rep, l);
		else mpz_add_ui( (mpz_ptr)&(res.gmp_rep), (mpz_srcptr)&gmp_rep, -l);
		return res;
	}

	// -- operator -
	Integer operator - (const int l, const Integer& n)
	{
		return -(n - (long)l);
	}
	Integer operator - (const unsigned int l, const Integer& n)
	{
		return -(n - (unsigned long)l);
	}
	Integer operator - (const long l, const Integer& n)
	{
		return -(n - l);
	}
	Integer operator - (const unsigned long l, const Integer& n)
	{
		return -(n - l);
	}
	Integer operator - (const Integer& n, const int l)
	{
		return n - (long)l;
	}
	Integer operator - (const Integer& n, const unsigned int l)
	{
		return n - (unsigned long)l;
	}

	Integer& operator -= (Integer& n, const int l)
	{
		return n -= (long)l;
	}
	Integer& operator -= (Integer& n, const unsigned int l)
	{
		return n -= (unsigned long)l;
	}

#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
	Integer operator - (const Integer& n, const long long l)
	{
		return n - (Integer)l;
	}
	Integer operator - (const Integer& n, const unsigned long long l)
	{
		return n - (Integer)l;
	}
	Integer operator - (const long long l, const Integer& n)
	{
		return n-l;
	}
	Integer operator - (const unsigned long long l, const Integer& n)
	{
		return n-l;
	}
	Integer& operator -= (Integer& n, const long long l)
	{
		return n -= (Integer)l;
	}
	Integer& operator -= (Integer& n, const unsigned long long l)
	{
		return n -= (Integer)l;
	}
#endif

}

#endif // __GIVARO_gmpxx_gmpxx_int_sub_C

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen
