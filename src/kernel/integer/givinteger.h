// =============================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier, J.-G. Dumas
// =============================================================

/*! @file givinteger.h
 * @ingroup integers
 * @brief Integer Domain class definition.
 */

#ifndef __GIVARO_integer_H
#define __GIVARO_integer_H

#include "givaro-config.h"
#include "gmp++/gmp++.h"
#include "givaro/givbasictype.h"
#include "givaro/givinit.h"
#include "givaro/giverror.h"
#include "givaro/givranditer.h"
#include "givaro/zring.h"
#include "givaro/ring-interface.h"
#include <string>

namespace Givaro {
	//------------------------------------ Class IntegerDom
	//! Integer Domain.
	class IntegerDom : public RingInterface<Integer> {
	public:
        using Self_t = IntegerDom;
		typedef Integer Rep;
		typedef Rep Element;


		IntegerDom() : one(1), mOne(-1), zero(0U) {}
		IntegerDom(const IntegerDom&) : one(1), mOne(-1), zero(0U) {}

		int operator==( const IntegerDom&) const
		{
			return 1;
		}
		int operator!=( const IntegerDom&) const
		{
			return 0;
		}

		// -- Constants:
		const Integer one;
		const Integer mOne;
		const Integer zero;

        Integer characteristic() const { return zero; }
        Integer& characteristic(Integer& p) const { return p = zero; }

		// -- assignement
		Rep& init  ( Rep& a ) const
		{
			return a;
		}
		Rep& init  ( Rep& a, const Rep& b) const
		{
			return a = b ;
		}
		Rep& read( Rep& a, const int64_t i) const
		{
			return a = Integer(i) ;
		}
		Rep& read( Rep& a, const uint64_t i) const
		{
			return a = Integer(i) ;
		}
		Rep& read( Rep& a, const int32_t i) const
		{
			return a = Integer(i) ;
		}
		Rep& read( Rep& a, const uint32_t i) const
		{
			return a = Integer(i) ;
		}

		Rep& convert( Rep& a, const Rep& b) const
		{
			return a = b ;
		}
		Rep& assign( Rep& a, const Rep& b) const
		{
			return a = b ;
		}

		// -- access
		const Rep& access(const Rep& a) const
		{
			return a;
		}

		template<class XXX> XXX& convert(XXX& x, const Rep& a) const
		{
			return Caster<XXX,Rep>(x,a);
		}

		// -- arithmetic operators
		Rep& mul( Rep& r, const Rep& a, const Rep& b ) const
		{
			return Integer::mul(r,a,b);
		}
		Rep& div( Rep& r, const Rep& a, const Rep& b ) const
		{
			return Integer::div(r,a,b);
		}
		Rep& mod( Rep& r, const Rep& a, const Rep& b ) const
		{
			return Integer::mod(r,a,b);
		}
		Rep& add( Rep& r, const Rep& a, const Rep& b ) const
		{
			return Integer::add(r,a,b);
		}
		Rep& sub( Rep& r, const Rep& a, const Rep& b ) const
		{
			return Integer::sub(r,a,b);
		}
		Rep& divmod( Rep& q, Rep& r, const Rep& a, const Rep& b ) const
		{
			return Integer::divmod(q,r,a,b);
		}
		Rep& divexact( Rep& q, const Rep& a, const Rep& b ) const
		{
			return Integer::divexact(q,a,b);
		}

		Rep& mulin( Rep& r, const Rep& a) const
		{
			return r *= a;
		}
		Rep& divin( Rep& r, const Rep& a) const
		{
			return r /= a;
		}
		Rep& modin( Rep& r, const Rep& a) const
		{
			return r %= a;
		}
		Rep& addin( Rep& r, const Rep& a) const
		{
			return r += a;
		}
		Rep& subin( Rep& r, const Rep& a) const
		{
			return r -= a;
		}

		Rep& axpy( Rep& r, const Rep& a, const Rep& b, const Rep& c ) const
		{
			return Integer::axpy(r,a,b,c);
		}
		Rep& maxpy( Rep& r, const Rep& a, const Rep& b, const Rep& c ) const
		{
			return Integer::maxpy(r,a,b,c);
		}
		Rep& maxpyin( Rep& r, const Rep& a, const Rep& b) const
		{
			return Integer::maxpyin(r,a,b);
		}
		Rep& axmy( Rep& r, const Rep& a, const Rep& b, const Rep& c ) const
		{
			return Integer::axmy(r,a,b,c);
		}
		Rep& axpyin( Rep& r, const Rep& a, const Rep& b ) const
		{
			return Integer::axpyin(r,a,b);
		}
		Rep& axmyin( Rep& r, const Rep& a, const Rep& b ) const
		{
			return Integer::axmyin(r,a,b);
		}

		// -- unary methods
		Rep& neg( Rep& r, const Rep& a ) const
		{
			return Integer::neg(r,a);
		}
		Rep& negin( Rep& r ) const
		{
			return Integer::negin(r);
		}

		// -- extended gcd  q = gcd(a,b) = u*a+v*b;
		Rep& gcd( Rep& g, Rep& u, Rep& v, const Rep& a, const Rep& b ) const
		{
			return ::Givaro::gcd(g, u, v, a, b);
		}
		Rep& gcd( Rep& g, const Rep& a, const Rep& b ) const
		{
			return ::Givaro::gcd(g, a, b);
		}
		Rep& gcdin( Rep& g, const Rep& a) const
		{
            Rep tmp(g);
			return ::Givaro::gcd(g, tmp, a);
		}
		Rep& lcm( Rep& l, const Rep& a, const Rep& b ) const
		{
			return ::Givaro::lcm(l, a, b);
		}
		Rep& lcmin( Rep& l, const Rep& a) const
		{ Rep tmp(l); return lcm(l, tmp, a);
		}

		Rep& inv(Rep& u, const Rep& a, const Rep& b) const
		{
			return ::Givaro::inv(u,a,b);
		}
		Rep& invin(Rep& u, const Rep& b) const
		{
			return ::Givaro::invin(u,b);
		}


		// - return n^l
		Rep& pow(Rep& r, const Rep& n, const int64_t l) const
		{
			return r = ::Givaro::pow(n, l);
		}
		Rep& pow(Rep& r, const Rep& n, const uint64_t l) const
		{
			return r = ::Givaro::pow(n, l);
		}
		Rep& pow(Rep& r, const Rep& n, const int32_t l) const
		{
			return r = ::Givaro::pow(n, (int64_t)l);
		}
		Rep& pow(Rep& r, const Rep& n, const uint32_t l) const
		{
			return r = ::Givaro::pow(n, (uint64_t)l);
		}

		// - return square root of n
		Rep& sqrt(Rep& s, const Rep& n) const
		{
			return ::Givaro::sqrt(s,n);
		}
		Rep& sqrt(Rep& s, Rep& r, const Rep& n) const
		{
			return ::Givaro::sqrtrem(s,n, r);
		}
		// - base p logarithm of a
		int64_t logp(const Rep& a, const Rep& p) const
		{
			return ::Givaro::logp(a,p);
		}

		// - return n^e % m
		Rep& powmod(Rep& r, const Rep& n, const int64_t e, const Rep& m) const
		{
			return r = ::Givaro::powmod(n, e, m);
		}
		Rep& powmod(Rep& r, const Rep& n, const Rep& e, const Rep& m) const
		{
			return r = ::Givaro::powmod(n, e, m);
		}

		// - Misc
		uint64_t length (const Rep& a) const
		{
			return ::Givaro::length(a);
		}
		int sign   (const Rep& a) const
		{
			return ::Givaro::sign(a);
		}
		bool isZero (const Rep& a) const
		{
			return ::Givaro::isZero(a);
		}
		bool isOne  (const Rep& a) const
		{
			return ::Givaro::isOne(a);
		}
		bool isMOne  (const Rep& a) const
		{
			return ::Givaro::isMOne(a);
		}
		/// isUnit
		inline  bool isUnit (const Rep& x) const
		{
			return ::Givaro::isOne(x) || ::Givaro::isMOne(x);
		}

		/** @brief isDivisor (a, b)
		 *  Test if b | a.
		 */
		inline  bool isDivisor (const Element& a, const Element& b) const
		{
			Element r; return ::Givaro::isZero(mod(r,a,b));
		}

        Element& abs(Element& x, const Element& a) const {
            return x=::Givaro::abs(a);
        }

        int32_t compare(const Rep& a, const Rep& b) const {
            return ::Givaro::compare(a,b);
        }
                
		bool areEqual (const Rep& a, const Rep& b) const
		{
			return ::Givaro::compare(a,b) ==0;
		}
		bool areNEqual(const Rep& a, const Rep& b) const
		{
			return ::Givaro::compare(a,b) !=0;
		}
		bool isgeq(const Rep& a, const Rep& b) const
		{
			return ::Givaro::compare(a,b) >= 0;
		}
		bool isleq(const Rep& a, const Rep& b) const
		{
			return ::Givaro::compare(a,b) <= 0;
		}
		bool isgeq(const int64_t b,const Rep& a ) const
		{
			return this->isgeq(Rep(b),a);
		}
		bool isleq(const int64_t b,const Rep& a ) const
		{
			return this->isleq(Rep(b),a);
		}
		bool isgeq(const Rep& a, const int64_t b) const
		{
			return this->isgeq(a,Rep(b));
		}
		bool isleq(const Rep& a, const int64_t b) const
		{
			return this->isleq(a,Rep(b));
		}
		bool isgt(const Rep& a, const Rep& b) const
		{
			return ::Givaro::compare(a,b) > 0;
		}
		bool islt(const Rep& a, const Rep& b) const
		{
			return ::Givaro::compare(a,b) < 0;
		}
		bool isgt(const int64_t b,const Rep& a ) const
		{
			return this->isgt(Rep(b),a);
		}
		bool islt(const int64_t b,const Rep& a ) const
		{
			return this->islt(Rep(b),a);
		}
		bool isgt(const Rep& a, const int64_t b) const
		{
			return this->isgt(a,Rep(b));
		}
		bool islt(const Rep& a, const int64_t b) const
		{
			return this->islt(a,Rep(b));
		}


#ifdef __GMP_PLUSPLUS__
		void seeding(unsigned long s = 0) const
		{ Integer::seeding(s) ;
		}
#endif

        typedef GeneralRingRandIter<Self_t> RandIter;
        typedef GeneralRingNonZeroRandIter<Self_t> NonZeroRandIter;


		template< class MyRandIter > Rep& random(MyRandIter&, Rep& r, long s = 1) const
		{
			return Integer::random(r,s);
		}
		template< class MyRandIter > Rep& random(MyRandIter&, Rep& r, const Rep& b) const
		{
			return Integer::random(r,b);
		}
		template< class MyRandIter > Rep& nonzerorandom(MyRandIter&, Rep& r, long s = 1) const
		{
			return Integer::nonzerorandom(r,s);
		}
		template< class MyRandIter > Rep& nonzerorandom (MyRandIter&,Rep& r, const Rep& b) const
		{
			return Integer::nonzerorandom(r,b);
		}

		// -- IO
		std::istream& read ( std::istream& i )
		{
			char ch;
			i >> std::ws >> ch;
			// JGD 22.03.03
			//    if (ch != 'I')
			//      GivError::throw_error(GivBadFormat("IntegerDom::read: bad signature domain"));
			return i;
		}
		std::ostream& write( std::ostream& o ) const
		{
			return o << 'I';
		}
		std::istream& read ( std::istream& i, Rep& n) const
		{
			return i >> n;
		}
		std::ostream& write( std::ostream& o, const Rep& n) const
		{
			return o << n;
		}
	};

} // Givaro

#endif //__GIVARO_integer_H
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
