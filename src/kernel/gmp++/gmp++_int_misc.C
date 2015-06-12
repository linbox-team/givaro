// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int_misc.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: M. Samama, T. Gautier
// $Id: gmp++_int_misc.C,v 1.16 2010-12-16 16:54:38 jgdumas Exp $
// ==========================================================================
// Description:
/** @file gmp++/gmp++_int_misc.C
 * miscing stuff.
 */

#ifndef __GIVARO_gmpxx_gmpxx_int_misc_C
#define __GIVARO_gmpxx_gmpxx_int_misc_C

#include <iostream>
#include <cmath>
#ifndef __GIVARO_INLINE_ALL
#include "gmp++/gmp++.h"
#endif
#ifndef __GIVARO_GMP_NO_CXX
#include <sstream>
#endif

namespace Givaro {
	//-------------------------------------------fact (uint64_t l)
	Integer fact ( uint64_t l)
	{
		Integer Res ;
		mpz_fac_ui( (mpz_ptr)&(Res.gmp_rep), l ) ;
		return Res ;
	}

	//-------------------------------------------square root
	Integer& sqrt(Integer& q, const Integer &a)
	{
		mpz_sqrt( (mpz_ptr)&(q.gmp_rep),
			  (mpz_srcptr)&(a.gmp_rep)) ;
		return q;
	}

	Integer& sqrtrem(Integer& q, const Integer &a, Integer& r)
	{
		mpz_sqrtrem( (mpz_ptr)&(q.gmp_rep),
			     (mpz_ptr)&(r.gmp_rep), (mpz_srcptr)&(a.gmp_rep)) ;
		return q;
	}

	Integer sqrt(const Integer &a)
	{
		Integer q;
		return sqrt(q,a);
	}

	Integer sqrtrem(const Integer &a, Integer& r)
	{
		Integer q;
		return sqrtrem(q,a,r);
	}

	bool root(Integer& q, const Integer &a, uint32_t n)
	{
		return (bool)mpz_root ((mpz_ptr)&(q.gmp_rep),
				       (mpz_srcptr)&(a.gmp_rep),
				       n);
	}

	void swap(Integer& a, Integer& b)
	{
		return mpz_swap( (mpz_ptr)&(a.gmp_rep), (mpz_ptr)&(b.gmp_rep));
	}


	// Natural logarithm of a
	// log(2) being close to 0.69314718055994531
	double naturallog(const Integer& a)
	{
		long int exp_;
		double d = mpz_get_d_2exp( &exp_, (mpz_srcptr)&(a.gmp_rep) );
		return (double)exp_*0.69314718055994531+log(d);
	}


	/*! Tests parity of an integer
	 * @param a integer
	 * @return 1 if odd, 0 if even
	 */
	bool isOdd(const Integer &a)
	{
		int32_t o = mpz_tstbit( (mpz_srcptr) &(a.gmp_rep), 0);
		return (o!=0); // or maybe should I write l==1 ^^
	}

	// base p logarithm of a
	int64_t logp(const Integer& a, const Integer& p)
	{
		std::list< Integer > pows;
		Integer puiss = p, sq;
		do {
			pows.push_back( puiss );
		} while ( (puiss *= puiss) <= a );
		puiss = pows.back(); pows.pop_back();
		int64_t res = (1 << pows.size());
		while (! pows.empty() ) {
			if ((sq = puiss * pows.back()) <= a) {
				puiss = sq;
				pows.pop_back();
				res += (1 << pows.size());
			} else
				pows.pop_back();
		}
		return res;
	}

	// approximation of the base 2 logarithm of a
	// 1/log(2) being close to 1.44269504088896341
	double logtwo(const Integer& a)
	{
		int64_t exp;
		double d = mpz_get_d_2exp( &exp, (mpz_srcptr)&(a.gmp_rep) );
		return (double)exp+log(d)*1.44269504088896341;
	}

	//------------------------------------------GMP isprime
	//     If this function returns 0, OP is definitely not prime.  If it
	//     returns 1, then OP is `probably' prime.  The probability of a
	//     false positive is (1/4)^r.  A reasonable value of r is 25.

	Integer& nextprime(Integer& r, const Integer &p)
	{
		mpz_nextprime ((mpz_ptr)&(r.gmp_rep), (mpz_srcptr)&(p.gmp_rep)) ;
		return r;
	}

	// Copied and adapted from mpz/nextprime.c
	Integer& prevprime(Integer& r, const Integer &p)
	{
		if (p < 3) return (r=2);
		if (isOdd(p))
			mpz_sub_ui ( (mpz_ptr)&(r.gmp_rep), (mpz_srcptr)&(p.gmp_rep), 2L );
		else
			mpz_sub_ui ( (mpz_ptr)&(r.gmp_rep), (mpz_srcptr)&(p.gmp_rep), 1 );

		while( !mpz_probab_prime_p ( (mpz_srcptr)&(p.gmp_rep), 10 ) )
			mpz_sub_ui ( (mpz_ptr)&(r.gmp_rep), (mpz_srcptr)&(p.gmp_rep), 2L );

		return r;
	}

	int32_t probab_prime(const Integer &p)
	{
		return mpz_probab_prime_p ((mpz_srcptr)&(p.gmp_rep),10) ;
	}

	int32_t probab_prime(const Integer &p, int32_t r)
	{
		return mpz_probab_prime_p ((mpz_srcptr)&(p.gmp_rep),r) ;
	}

	// ==========================================================================
	// Computes and returns the Jacobi and Legendre symbols (u/v) of the integers u and v.
	// The algorithm used is Gmp's.
	int32_t jacobi(const Integer& u, const Integer& v)
	{
		return mpz_jacobi ((mpz_srcptr)&(u.gmp_rep),(mpz_srcptr)&(v.gmp_rep)) ;
	}

	int32_t legendre(const Integer& u, const Integer& v)
	{
		return mpz_legendre ((mpz_srcptr)&(u.gmp_rep),(mpz_srcptr)&(v.gmp_rep)) ;
	}



	//--------------------------------------------Integer::operator <<   // shift left
	Integer Integer::operator << (int32_t l) const
	{
		return this->operator<<( (uint64_t)l );
	}
	Integer Integer::operator << (uint32_t l) const
	{
		return this->operator<<( (uint64_t)l );
	}
	Integer Integer::operator << (int64_t l) const
	{
		return this->operator<<( (uint64_t)l );
	}

	Integer Integer::operator << (uint64_t l) const
	{
		Integer tmp;
		mpz_mul_2exp((mpz_ptr)&(tmp.gmp_rep), (mpz_srcptr)&(gmp_rep), l );
		return tmp;
	}


	//--------------------------------------------Integer::operator >>   // shift right
	Integer Integer::operator >> (int32_t l) const
	{
		return this->operator>>( (uint64_t)l );
	}

	Integer Integer::operator >> (int64_t l) const
	{
		return this->operator>>( (uint64_t)l );
	}

	Integer Integer::operator >> (uint32_t l) const
	{
		return this->operator>>( (uint64_t)l );
	}

	Integer Integer::operator >> (uint64_t l) const
	{
		Integer tmp;
		mpz_tdiv_q_2exp( (mpz_ptr)&(tmp.gmp_rep), (mpz_srcptr)&(gmp_rep), l );
		return tmp;
	}

	//--------------------------------------------Integer::operator <<=   // shift left
	Integer& Integer::operator <<= (int32_t l)
	{
		return this->operator<<= ( (uint64_t)l );
	}
	Integer& Integer::operator <<=  (uint32_t l)
	{
		return this->operator<<= ( (uint64_t)l );
	}
	Integer& Integer::operator <<= (int64_t l)
	{
		return this->operator<<= ( (uint64_t)l );
	}

	Integer& Integer::operator <<= (uint64_t l)
	{
		mpz_mul_2exp((mpz_ptr)&(gmp_rep), (mpz_srcptr)&(gmp_rep), l );
		return *this;
	}


	//--------------------------------------------Integer::operator >>=   // shift right
	Integer& Integer::operator >>= (int32_t l)
	{
		return this->operator>>= ( (uint64_t)l );
	}
	Integer& Integer::operator >>= (int64_t l)
	{
		return this->operator>>= ( (uint64_t)l );
	}
	Integer& Integer::operator >>= (uint32_t l)
	{
		return this->operator>>= ( (uint64_t)l );
	}

	Integer& Integer::operator >>= (uint64_t l)
	{
		mpz_tdiv_q_2exp( (mpz_ptr)&(gmp_rep), (mpz_srcptr)&(gmp_rep), l );
		return *this;
	}

	//------------------------------------------- Bit logic
	Integer Integer::operator^ (const Integer& a) const
	{   // XOR
		Integer res(*this);
		return res ^= a;
	}
	Integer Integer::operator| (const Integer& a) const
	{   // OR
		Integer res(*this);
		return res |= a;
	}
	Integer Integer::operator& (const Integer& a) const
	{   // AND
		Integer res(*this);
		return res &= a;
	}
	Integer Integer::operator^ (const uint64_t & a) const
	{   // XOR
		Integer res(*this);
		return res ^= a;
	}
	Integer Integer::operator| (const uint64_t & a) const
	{   // OR
		Integer res(*this);
		return res |= a;
	}
	uint64_t Integer::operator& (const uint64_t & a) const
	{   // AND
		return mpz_get_ui((mpz_srcptr)&(gmp_rep)) & a;
	}
	Integer Integer::operator^ (const uint32_t& a) const
	{   // XOR
		Integer res(*this);
		return res ^= a;
	}
	Integer Integer::operator| (const uint32_t& a) const
	{   // OR
		Integer res(*this);
		return res |= a;
	}
	uint32_t Integer::operator& (const uint32_t& a) const
	{   // AND
		return (uint32_t) (mpz_get_ui((mpz_srcptr)&(gmp_rep)) & (uint64_t)a );
	}
	Integer Integer::operator~ () const
	{   // 1 complement
		Integer res;
		mpz_com( (mpz_ptr)&(res.gmp_rep), (mpz_srcptr)&(gmp_rep));
		return res;
	}
	Integer& Integer::operator^= (const Integer& a)
	{   // XOR
		mpz_xor( (mpz_ptr)&(gmp_rep), (mpz_ptr)&(gmp_rep), (mpz_srcptr)&(a.gmp_rep));
		return *this;
	}
	Integer& Integer::operator|= (const Integer& a)
	{   // OR
		mpz_ior( (mpz_ptr)&(gmp_rep), (mpz_ptr)&(gmp_rep), (mpz_srcptr)&(a.gmp_rep));
		return *this;
	}
	Integer& Integer::operator&= (const Integer& a)
	{   // AND
		mpz_and( (mpz_ptr)&(gmp_rep), (mpz_ptr)&(gmp_rep), (mpz_srcptr)&(a.gmp_rep));
		return *this;
	}

	Integer& Integer::operator^= (const uint64_t & a)
	{   // XOR
        Integer au(a);
		mpz_xor( (mpz_ptr)&(gmp_rep), (mpz_ptr)&(gmp_rep), (mpz_srcptr)&(au.gmp_rep));
		return *this;
	}
	Integer& Integer::operator|= (const uint64_t & a)
	{   // OR
        Integer au(a);
		mpz_ior( (mpz_ptr)&(gmp_rep), (mpz_ptr)&(gmp_rep), (mpz_srcptr)&(au.gmp_rep));
		return *this;
	}
	Integer& Integer::operator&= (const uint64_t & a)
	{   // AND
        Integer au(a);
		mpz_and( (mpz_ptr)&(gmp_rep), (mpz_ptr)&(gmp_rep), (mpz_srcptr)&(au.gmp_rep));
		return *this;
	}

	Integer& Integer::operator^= (const uint32_t& a)
	{   // XOR
        Integer au(a);
		mpz_xor( (mpz_ptr)&(gmp_rep), (mpz_ptr)&(gmp_rep), (mpz_srcptr)&(au.gmp_rep));
		return *this;
	}
	Integer& Integer::operator|= (const uint32_t& a)
	{   // OR
        Integer au(a);
		mpz_ior( (mpz_ptr)&(gmp_rep), (mpz_ptr)&(gmp_rep), (mpz_srcptr)&(au.gmp_rep));
		return *this;
	}
	Integer& Integer::operator&= (const uint32_t& a)
	{   // AND
        Integer au(a);
		mpz_and( (mpz_ptr)&(gmp_rep), (mpz_ptr)&(gmp_rep), (mpz_srcptr)&(au.gmp_rep));
		return *this;
	}


	//------------------------------------------- convert method
	//------------------------------------------- casting method
	Integer::operator int() const
	{
		return int32_t (mpz_get_si ( (mpz_srcptr)&gmp_rep));
	}
	Integer::operator uint32_t() const
	{
		return (uint32_t) mpz_get_ui ( (mpz_srcptr)&gmp_rep);
	}
	Integer::operator int64_t() const
	{
		return mpz_get_si ( (mpz_srcptr)&gmp_rep);
	}
	Integer::operator uint64_t() const
	{
		return mpz_get_ui ( (mpz_srcptr)&gmp_rep);
	}
	Integer::operator double() const
	{
		return mpz_get_d ( (mpz_srcptr)&gmp_rep);
	}
	Integer::operator float() const
	{
		return (float)mpz_get_d ( (mpz_srcptr)&gmp_rep);
	}

	Integer::operator std::string () const
	{
#ifdef __GIVARO_GMP_NO_CXX
		std::string s;
		uint64_t strSize = mpz_sizeinbase((mpz_srcptr)&(gmp_rep), 10) + 2;
		char *str = new char[strSize + 2];
		mpz_get_str(str, 10, (mpz_srcptr)&(gmp_rep));
		s = std::string(str);
		delete [] str ;
		// return ??
#else
		std::ostringstream o ;
		print(o);
		return o.str();
#endif
	}

	Integer::operator Integer::vect_t () const
	{
		size_t s = mpz_size( (mpz_srcptr)&(gmp_rep) );
		std::vector<mp_limb_t> v(s);
		std::vector<mp_limb_t>::iterator vi = v.begin();
		for(mp_size_t i = 0;vi != v.end();++vi, ++i)
			*vi = mpz_getlimbn( (mpz_srcptr)& (gmp_rep) ,i);
		return v;
	}

	uint64_t length(const Integer& a)
	{
		//! @bug JGD 23.04.2012: shouldn't it be "mp_limb_t" instead of "uint64_t"?
		return mpz_size( (mpz_srcptr)&(a.gmp_rep) ) * sizeof(uint64_t);
	}

	Integer abs(const Integer &n)
	{
		if (sign(n) >= 0)
			return n;
		return -n;
	}

	size_t Integer::size() const
	{
		return  mpz_size( (mpz_srcptr)&gmp_rep ) ;
	}

	size_t Integer::size_in_base(int32_t BASE) const
	{
		return  mpz_sizeinbase ((mpz_srcptr)&gmp_rep, BASE);
	}

	size_t Integer::bitsize() const
	{
		return  mpz_sizeinbase ((mpz_srcptr)&gmp_rep, 2);
	}

	uint64_t Integer::operator[](size_t i) const
	{
		if ( mpz_size( (mpz_srcptr)&gmp_rep ) > i)
			return mpz_getlimbn( (mpz_srcptr)&gmp_rep, i);
		else
			return 0;
	}

}

#endif // __GIVARO_gmpxx_gmpxx_int_misc_C

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen
