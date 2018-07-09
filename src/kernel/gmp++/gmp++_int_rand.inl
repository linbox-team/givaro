// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int_add.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: B. Boyer
// $Id: gmp++_int_random.C,v 1.5 2011-09-16 12:09:37 briceboyer Exp $
// ==========================================================================

// #include "gmp++/gmp++.h"
/** @file gmp++/gmp++_int_rand.inl
 * randing stuff.
 */

#ifndef __GIVARO_gmpxx_gmpxx_int_rand_INL
#define __GIVARO_gmpxx_gmpxx_int_rand_INL

#include <givaro/givtimer.h>
#include <givaro/givrandom.h>
#include <givaro/udl.h>
#include <assert.h>

namespace Givaro {
	//-----------------------------------------------------
	//----------------------- Random integers -------------
	//-----------------------------------------------------

	template <class RNG>
	bool Integer::RandBool(RNG& rs)
	{
		return (rs() & 1u);
	}


	/* ****************************** */
	/*  random number smaller than m  */
	/* ****************************** */

	//! returns a random integer \p r in the intervall <code>[[x, m-1]]</code>
	//! where x = 0 or -(m-1) according to \p ALWAYSPOSITIVE
	//! @bug m \b has to be an integer here.
	//@{
	template<bool ALWAYSPOSITIVE, class RNG>
	Integer& Integer::random_lessthan (Integer& r, const Integer & m, RNG& rs)
	{
		/*
		mpz_set( (mpz_ptr) &(r.gmp_rep) ,
			 ( (mpz_class)rs.getIntRand().get_z_range((mpz_class) (mpz_srcptr) &(m.gmp_rep)) ).get_mpz_t() );
		*/
		rs.getIntRand().get_z_range(r, m);
		if(!ALWAYSPOSITIVE) if (Integer::RandBool(rs)) Integer::negin(r);
		return r;
	}

	template <class RNG>
	Integer& Integer::random_lessthan (Integer& r, const Integer & m, RNG& rs)
	{
		return random_lessthan<true>(r,m,rs);
	}
	//@}

	/* ******************************** */
	/*  random number smaller than 2^m  */
	/* ******************************** */
	//! returns a random integer \p r in the intervall <code>[[x, 2^m-1]]</code>
	//! where x = 0 or -(2^m-1) according to \p ALWAYSPOSITIVE
	//! returns a random integer \p r of at most \p m bits
	//@{

	template<bool ALWAYSPOSITIVE, class RNG>
	Integer& Integer::random_lessthan_2exp (Integer& r, const uint64_t & m, RNG& rs)
	{
		//mpz_set( (mpz_ptr) &(r.gmp_rep) , ((mpz_class)rs.getIntRand().get_z_bits(m)).get_mpz_t() );
		rs.getIntRand().get_z_bits(r, m);
		if(!ALWAYSPOSITIVE) {
			if (Integer::RandBool(rs)) Integer::negin(r);
		}
		return r;
	}

	template<bool ALWAYSPOSITIVE, class RNG>
	Integer Integer::random_lessthan_2exp (const uint64_t & m, RNG& rs)
	{
		Integer r ;
		random_lessthan_2exp<ALWAYSPOSITIVE>(r,m,rs);
		return r;
	}

	/* synonyms CAREFULL: when m is integer, meaning is different*/
	template<bool ALWAYSPOSITIVE, class RNG>
	Integer& Integer::random_lessthan (Integer& r, const uint64_t & m, RNG& rs)
	{
		return Integer::random_lessthan_2exp<ALWAYSPOSITIVE>(r,m,rs);
	}
	template <class RNG>
	Integer& Integer::random_lessthan (Integer& r, const uint64_t & m, RNG& rs)
	{
		return random_lessthan<true>(r,m,rs);
	}

	template<bool ALWAYSPOSITIVE,class T, class RNG>
	Integer Integer::random_lessthan (const T & m, RNG& rs)
	{
		Integer res ;
		return Integer::random_lessthan<ALWAYSPOSITIVE>(res,(typename Signed_Trait<T>::unsigned_type)m,rs);
	}

	template <class RNG>
	Integer& Integer::random_lessthan_2exp (Integer& r, const uint64_t & m, RNG& rs)
	{
		return random_lessthan_2exp<true>(r,m,rs);
	}

	template <class RNG>
	Integer Integer::random_lessthan_2exp (const uint64_t & m, RNG& rs)
	{
		return random_lessthan_2exp<true>(m,rs);
	}

	template<class T, class RNG>
	Integer Integer::random_lessthan (const T & m, RNG& rs)
	{
		return random_lessthan<true>(m,rs);
	}
	//@}


	/* ********************************* */
	/*  random number of same size as s  */
	/* ********************************* */

	//! returns a reference to a random number \p r of the size of \p s, exactly.
	template<bool ALWAYSPOSITIVE, class RNG>
	Integer& Integer::random_exact (Integer& r, const Integer & s, RNG& rs)
	{
		size_t t = s.bitsize() ;
		Integer::random_exact_2exp<ALWAYSPOSITIVE>(r,t,rs);
		return r;
	}
	template <class RNG>
	Integer& Integer::random_exact (Integer& r, const uint64_t & m, RNG& rs)
	{
		return Integer::random_exact<true>(r,m,rs);
	}
	template <class RNG>
	Integer& Integer::random_exact (Integer& r, const Integer & s, RNG& rs)
	{
		return Integer::random_exact<true>(r,s,rs);
	}

	template<class T, class RNG>
	typename std::enable_if<IsGivRand<RNG>::value,Integer>::type
	Integer::random_exact (const T & s, RNG& rs)
	{
		return Integer::random_exact<true>(s,rs) ;
	}



	/* ************************* */
	/*  random number of size m  */
	/* ************************* */

	//! returns a reference to a random number \p r of the size \p m bits, exactly.
	template<bool ALWAYSPOSITIVE, class RNG>
	Integer& Integer::random_exact_2exp (Integer& r, const uint64_t & m, RNG& rs)
	{
		if (m) random_lessthan_2exp<true>(r,m-1_ui64,rs);
		mpz_setbit( (mpz_ptr) &(r.gmp_rep) , m-1_ui64);
		if(!ALWAYSPOSITIVE) if (Integer::RandBool(rs)) Integer::negin(r);
		return r;
	}

	template <class RNG>
	Integer& Integer::random_exact_2exp (Integer& r, const uint64_t & m, RNG& rs)
	{
		return Integer::random_exact_2exp<true>(r,m,rs);
	}
	// synonym
	template<bool ALWAYSPOSITIVE, class RNG>
	Integer& Integer::random_exact (Integer& r, const uint64_t & m, RNG& rs)
	{
		return Integer::random_exact_2exp<ALWAYSPOSITIVE>(r,m,rs) ;
	}

	template<bool ALWAYSPOSITIVE,class T, class RNG>
	typename std::enable_if<IsGivRand<RNG>::value,Integer>::type
	Integer::random_exact (const T & s, RNG& rs)
	{
		Integer res ;
		return random_exact<ALWAYSPOSITIVE>(res,s,rs);
	}

	/* **************************** */
	/*  random number in [[m,M-1]]  */
	/* **************************** */

	template <class RNG>
	Integer& Integer::random_between (Integer& r, const Integer& m, const Integer&M, RNG& rs)
	{
		assert(M > m);
		random_lessthan(r,Integer(M-m),rs);
		r += m ;
		return (r);
	}

	template <class RNG>
	Integer Integer::random_between (const Integer& m, const Integer &M, RNG& rs)
	{
		Integer r ;
		return random_between(r,m,M,rs);
	}


	template<class R, class RNG>
	typename std::enable_if<IsGivRand<RNG>::value,Integer>::type
	Integer::random_between (const R & m, const R & M, RNG& rs)
	{
		return Integer::random_between(static_cast<uint64_t>(m),
					       static_cast<uint64_t>(M),
					       rs);
	}
	template<class R, class RNG>
	typename std::enable_if<!IsGivRand<R>::value,Integer>::type&
	Integer::random_between (Integer &r, const R & m, const R & M, RNG& rs)
	{
		return Integer::random_between(r,static_cast<uint64_t>(m),
					       static_cast<uint64_t>(M),
					       rs);
	}


	/* ******************************** */
	/*  random number in [[2^m,2^M-1]]  */
	/* ******************************** */
	// todo : template<bool ALWAYSPOSITIVE, bool V, class RNG>
	template <class RNG>
	Integer& Integer::random_between_2exp (Integer& r, const uint64_t& m, const uint64_t &M, RNG& rs)
	{
		assert(M > m);
		r = nonzerorandom((uint64_t)M-m,rs);
		Integer r1 = random_lessthan_2exp(m,rs);
		r <<= m ;
		r+= r1 ;
		return (r);
	}

	template <class RNG>
	Integer Integer::random_between_2exp (const uint64_t & m, const uint64_t &M, RNG& rs)
	{
		Integer r ;
		return random_between_2exp(r,m,M,rs);
	}

	// synonym.
	template <class RNG>
	Integer Integer::random_between (const uint64_t & m, const uint64_t &M, RNG& rs)
	{
		return random_between_2exp(m,M,rs) ;
	}


	template <class RNG>
	Integer& Integer::random_between (Integer& r, const uint64_t& m, const uint64_t &M, RNG& rs)
	{
		return random_between_2exp(r,m,M,rs);
	}
	/* **************/
	/*  short hand  */
	/* **************/

	//! returns a random integer less than...
	template<bool ALWAYSPOSITIVE,class T, class RNG>
	typename std::enable_if<!IsGivRand<T>::value,Integer>::type&
	Integer::random (Integer& r, const T & m, RNG& rs)
	{
		return Integer::random_lessthan<ALWAYSPOSITIVE>(r, (typename Signed_Trait<T>::unsigned_type) m, rs) ;
	}

	//! returns a random integer less than...
	template<bool ALWAYSPOSITIVE,class T, class RNG>
	typename std::enable_if<IsGivRand<RNG>::value,Integer>::type
	Integer::random(const T & sz, RNG& rs)
	{
		return Integer::random_lessthan<ALWAYSPOSITIVE,T>(sz,rs);
	}

	template <class RNG>
	typename std::enable_if<IsGivRand<RNG>::value,Integer>::type
	Integer::random(RNG& rs)
	{
		return Integer::random(sizeof(mp_limb_t)*8,rs) ;
	}
	template<bool ALWAYSPOSITIVE, class RNG>
	typename std::enable_if<IsGivRand<RNG>::value,Integer>::type
	Integer::random(RNG& rs)
	{
		Integer rez = Integer::random(sizeof(mp_limb_t)*8,rs) ;
		if (!ALWAYSPOSITIVE) if (Integer::RandBool(rs)) negin(rez);
		return rez;
	}

	template<class T, class RNG>
	typename std::enable_if<!IsGivRand<T>::value,Integer>::type&
	Integer::random (Integer& r, const T & m, RNG& rs)
	{
		return Integer::random<true,T>(r,m,rs);
	}

	template<class T, class RNG>
	typename std::enable_if<IsGivRand<RNG>::value,Integer>::type
	Integer::random(const T & sz, RNG& rs)
	{
		return Integer::random<true>(sz,rs);
	}

	/* *******************/
	/*  Non Zero random  */
	/* *******************/

	template<bool ALWAYSPOSITIVE, class T, class RNG>
	typename std::enable_if<IsGivRand<RNG>::value,Integer>::type
	Integer::nonzerorandom(const T & sz, RNG& rs)
	{
		Integer r;
		while(isZero(Integer::random<ALWAYSPOSITIVE,T>(r, sz, rs) )) {} ;
		return r;
	}

	// BB: It's also 1+random(sz-1)...

	template<bool ALWAYSPOSITIVE, class T, class RNG>
	typename std::enable_if<!IsGivRand<T>::value,Integer>::type&
	Integer::nonzerorandom (Integer& r, const T& size, RNG& rs)
	{
		while (isZero(Integer::random<ALWAYSPOSITIVE,T>(r,size,rs))) {} ;
		return r;
	}

	template<class T, class RNG>
	typename std::enable_if<IsGivRand<RNG>::value,Integer>::type
	Integer::nonzerorandom(const T & sz, RNG& rs)
	{
		return  Integer::nonzerorandom<true>(sz,rs);
	}
	template<class T, class RNG>
	typename std::enable_if<!IsGivRand<T>::value,Integer>::type&
	Integer::nonzerorandom (Integer& r, const T& size, RNG& rs)
	{
		return  Integer::nonzerorandom<true>(r,size,rs);
	}
	template <class RNG>
	Integer  Integer::nonzerorandom(RNG& rs)
	{
		Integer rez = Integer::nonzerorandom(sizeof(mp_limb_t)*8,rs) ;
		// if (!ALWAYSPOSITIVE) if (Integer::RandBool()) negin(rez);
		return rez;
	}

	// complete definitions of some functions in givrandom.h
	inline void GivIntRand::seed(const Integer& s)
	{
		gmp_randseed(_state, s.get_rep());
	}

	inline void GivIntRand::get_z_range(Integer& result, const Integer& m)
	{
		mpz_urandomm(result.get_rep(), _state, m.get_rep());
	}

	// result <$- {0,...,2^bits-1}
	inline void GivIntRand::get_z_bits(Integer& result, uint64_t bits)
	{
		mpz_urandomb(result.get_rep(), _state, bits);
	}

	template <class EngineType, class IntRNGType>
	void GivRandomGeneric<EngineType,IntRNGType>::seed(const Integer& s)
	{
		if (!s) seed(0);
		_engine.seed(s[0]);
		_intrng.seed(s);
	}

}

#endif // __GIVARO_gmpxx_gmpxx_int_rand_INL

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen
