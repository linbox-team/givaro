/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* Copyright (C) 2014 LinBox
* Written by BB <brice.boyer@lip6.fr>
*
* ========LICENCE========
* This file is part of the library FFLAS-FFPACK.
*
* FFLAS-FFPACK is free software: you can redistribute it and/or modify
* it under the terms of the  GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
* This library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with this library; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
* ========LICENCE========
*
*/

/*! @file field/modular-uint64.h
* @ingroup field
* @brief  representation of <code>Z/mZ</code> over \c uint64_t .
*/

#ifndef __GIVARO_modular_uint16_H
#define __GIVARO_modular_uint16_H

#include "givaro/givranditer.h"
#include "givaro/ring-interface.h"
#include "givaro/modular-general.h"

namespace Givaro {

/** \brief Specialization of Modular to uint64_t element type with efficient dot product.
 *
 * Efficient element operations for dot product, mul, axpy, by using floating point
 * inverse of modulus (borrowed from NTL) and some use of non-normalized intermediate values.
 *
 * For some uses this is the most efficient field for primes in the range from half word
 * to 2^30.
 *
 * Requires: Modulus < 2^30.
 * Intended use: 2^15 < prime modulus < 2^30.
 * \ingroup field
 */

template <typename COMP>
class Modular<uint16_t, COMP> : public RingInterface<uint16_t>
{
public:
	// ----- Exported Types and constantes
	using Self_t = Modular<uint16_t, COMP>;
	using Compute_t = typename std::make_unsigned<COMP>::type;
	using Residu_t = uint16_t;
	
	enum { size_rep = sizeof(Residu_t) };

	// ----- Constantes
	const Element zero;
	const Element one;
	const Element mOne;
	
	// ----- Constructors
	Modular()
	: zero(Element(0)), one(Element(1)), mOne(Element(0)), _p(Element(0)) {}

	Modular(Residu_t p)
	: zero(Element(0)), one(Element(1)), mOne(Element(p-1)), _p(p)
	{
	    assert(_p >= getMinModulus());
	    assert(_p <= getMaxModulus());
	}

	Modular(const Self_t& F)
	: zero(F.zero), one(F.one), mOne(F.mOne), _p(F._p) {}

	// ----- Accessors
	inline Element minElement() const final { return zero; }
	inline Element maxElement() const final { return mOne; }

	// ----- Access to the modulus
	inline Residu_t residu() const { return _p; }
	inline Residu_t size() const { return _p; }
	inline Residu_t characteristic() const { return _p; }
	inline Integer& characteristic(Integer& p) const { return p = _p; }
	inline Residu_t cardinality() const { return _p; }
	inline Integer& cardinality(Integer& p) const { return p = _p; }
	static inline Residu_t getMaxModulus();
	static inline Residu_t getMinModulus() { return 2; }

	// ----- Checkers
	inline bool isZero(const Element& a) const final { return a == zero; }
	inline bool isOne (const Element& a) const final { return a == one; }
	inline bool isMOne(const Element& a) const final { return a == mOne; }
	inline bool areEqual(const Element& a, const Element& b) const final { return a == b; }
	inline size_t length(const Element a) const { return size_rep; }
	
	// ----- Ring-wise operators
	bool operator==(const Self_t& F) const { return _p == F._p; }
	bool operator!=(const Self_t& F) const { return _p != F._p; }
	Self_t& operator=(const Self_t& F)
	{
		F.assign(const_cast<Element&>(one),  F.one);
		F.assign(const_cast<Element&>(zero), F.zero);
		F.assign(const_cast<Element&>(mOne), F.mOne);
		_p = F._p;
		return *this;
	}
	
	// ----- Initialisation
	Element& init (Element &x) const;
	Element& init (Element &x, const int32_t& y ) const;
	Element& init (Element &x, const int64_t& y ) const;
	Element& init (Element &x, const uint32_t& y ) const;
	Element& init (Element &x, const uint64_t& y ) const;
	Element& init (Element &x, const double& y) const;
	Element& init (Element& r, const Integer& a) const;
	template<class XXX> Element& init(Element& x, const XXX& y) const;

	Element &assign (Element &x, const Element &y) const;
	
	// ----- Convert
	Integer& convert(Integer& i, const Element a) const { unsigned long ur; return i = (Integer)convert(ur, a);	}
	template<typename XXX> XXX& convert(XXX& r, const Element a ) const { return r = static_cast<XXX>(a) ;}

	inline Element& reduce (Element& x, const Element& y) const { return x = y % _p; }
	inline Element& reduce (Element& x) const { return x %= _p; }

	// ----- Classic arithmetic
	Element& mul(Element& r, const Element& a, const Element& b) const final;
	Element& div(Element& r, const Element& a, const Element& b) const final;
	Element& add(Element& r, const Element& a, const Element& b) const final;
	Element& sub(Element& r, const Element& a, const Element& b) const final;
	Element& neg(Element& r, const Element& a) const final;
	Element& inv(Element& r, const Element& a) const final;

	Element& mulin(Element& r, const Element& a) const final;
	Element& divin(Element& r, const Element& a) const final;
	Element& addin(Element& r, const Element& a) const final;
	Element& subin(Element& r, const Element& a) const final;
	Element& negin(Element& r) const final;
	Element& invin(Element& r) const final;
	
	// -- axpy:   r <- a * x + y
	// -- axpyin: r <- a * x + r
	Element& axpy  (Element& r, const Element& a, const Element& x, const Element& y) const final;
	Element& axpyin(Element& r, const Element& a, const Element& x) const final;

	// -- axmy:   r <- a * x - y
	// -- axmyin: r <- a * x - r
	Element& axmy  (Element& r, const Element& a, const Element& x, const Element& y) const final;
	Element& axmyin(Element& r, const Element& a, const Element& x) const final;

	// -- maxpy:   r <- y - a * x
	// -- maxpyin: r <- r - a * x
	Element& maxpy  (Element& r, const Element& a, const Element& x, const Element& y) const final;
	Element& maxpyin(Element& r, const Element& a, const Element& x) const final;

	// ----- Random generators
	typedef ModularRandIter<Self_t> RandIter;
	typedef GeneralRingNonZeroRandIter<Self_t> NonZeroRandIter;
    template< class Random > Element& random(const Random& g, Element& r) const { return init(r, g()); }
    template< class Random > Element& nonzerorandom(const Random& g, Element& a) const
    	{ while (isZero(init(a, g())));
    	  return a; }

	// --- IO methods
	std::istream& read (std::istream& s);
	std::ostream& write(std::ostream& s) const;
	std::istream& read (std::istream& s, Element& a) const;
	std::ostream& write(std::ostream& s, const Element& a) const;

protected:
	// -- data representation of the domain:
	Residu_t _p;
};

}

#include "givaro/modular-uint16.inl"

#endif //__GIVARO_modular_uint16_H

