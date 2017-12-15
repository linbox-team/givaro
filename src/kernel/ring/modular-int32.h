/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givzpz32std.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givzpz32std.h,v 1.22 2011-02-04 14:11:46 jgdumas Exp $
// ==========================================================================
//
//  Modified by Pascal Giorgi on 2002/02/13  (pascal.giorgi@ens-lyon.fr)
//  Modified by Alexis Breust on 2014/12/16  (alexis.breust@imag.fr)
//  Modified by Romain Lebreton on 2016/06/10 (romain.lebreton@lirmm.fr)

/*! @file ring/modular-int32.h
 * @ingroup ring
 * @brief  representation of <code>Z/mZ</code> over \c int32_t .
 */

#ifndef __GIVARO_modular_int32_H
#define __GIVARO_modular_int32_H

#include "givaro/givinteger.h"
#include "givaro/givcaster.h"
#include "givaro/givranditer.h"
#include "givaro/ring-interface.h"
#include "givaro/modular-general.h"

namespace Givaro {

	/*! @brief This class implement the standard arithmetic with Modulo Elements.
	 * - The representation of an integer a in Zpz is the value a % p
	 * - m max is 46341
	 * - p max is 46337
	 * .
	 */
	template<typename COMP>
	class Modular<int32_t, COMP> : public FiniteFieldInterface<int32_t>
	{
	public:

		// ----- Exported Types and constantes
		using Self_t = Modular<int32_t, COMP>;
		using Residu_t = uint32_t;
		using Compute_t = typename std::make_unsigned<COMP>::type;
		enum { size_rep = sizeof(Residu_t) };

		// ----- Representation of vector of the Element
		typedef Element* Array;
		typedef const Element* constArray;

		// ----- Constantes
		const Element zero;
		const Element one;
		const Element mOne;

		// ----- Constructors
		Modular()
			: zero(static_cast<Element>(0))
			, one(static_cast<Element>(1))
			, mOne(static_cast<Element>(-1))
			, _p(static_cast<Residu_t>(0))
			, _bitsizep(0) {}

		Modular(const Residu_t p)
			: zero(static_cast<Element>(0))
			, one(static_cast<Element>(1))
			, mOne(static_cast<Element>(p-1))
			, _p(static_cast<Residu_t>(p))
			, _bitsizep(0)
		{
			assert(_p >= minCardinality());
			assert(_p <= maxCardinality());
			Residu_t __p = _p;
			while (__p != 0) {
				_bitsizep++;
				__p >>= 1;
			}
		}

		Modular(const Self_t& F)
			: zero(F.zero), one(F.one), mOne(F.mOne), _p(F._p), _bitsizep(F._bitsizep) {}

		// ----- Accessors
		inline Element minElement() const override { return zero; }
		inline Element maxElement() const override { return mOne; }

		// ----- Access to the modulus
		inline Residu_t residu() const { return _p; }
		inline Residu_t size() const { return _p; }
		inline Residu_t characteristic() const { return _p; }
		inline Residu_t cardinality() const { return _p; }
		template<class T> inline T& characteristic(T& p) const { return p = _p; }
		template<class T> inline T& cardinality(T& p) const { return p = _p; }

		static inline Residu_t maxCardinality();
		static inline Residu_t minCardinality() { return 2; }

		// ----- Checkers
		inline bool isZero(const Element& a) const override { return a == zero; }
		inline bool isOne (const Element& a) const override { return a == one; }
		inline bool isMOne(const Element& a) const override { return a == mOne; }
        inline bool isUnit(const Element& a) const override;
		inline bool areEqual(const Element& a, const Element& b) const override { return a == b; }
		inline size_t length(const Element a) const { return size_rep; }

		// ----- Ring-wise operators
		inline bool operator==(const Self_t& F) const { return _p == F._p; }
		inline bool operator!=(const Self_t& F) const { return _p != F._p; }
		inline Self_t& operator=(const Self_t& F)
		{
			F.assign(const_cast<Element&>(one),  F.one);
			F.assign(const_cast<Element&>(zero), F.zero);
			F.assign(const_cast<Element&>(mOne), F.mOne);
			_p = F._p;
			_bitsizep = F._bitsizep;
			return *this;
		}

		// ----- Initialisation
		Element& init (Element& x) const;
		Element& init (Element& x, const float y) const;
		Element& init (Element& x, const double y) const;
		Element& init (Element& x, const int64_t y) const;
		Element& init (Element& x, const uint64_t y) const;
		Element& init (Element& x, const Integer& y) const;
		template<typename T> Element& init(Element& r, const T& a) const
		{ return reduce(Caster<Element,T>(r,a)); }

		Element& assign (Element& x, const Element& y) const;

		// ----- Convert and reduce
		template<typename T> T& convert(T& r, const Element& a) const
		{ return Caster<T,Element>(r,a); }

		Element& reduce (Element& x, const Element& y) const;
		Element& reduce (Element& x) const;

		// ----- Classic arithmetic
		Element& mul(Element& r, const Element& a, const Element& b) const override;
		Element& div(Element& r, const Element& a, const Element& b) const override;
		Element& add(Element& r, const Element& a, const Element& b) const override;
		Element& sub(Element& r, const Element& a, const Element& b) const override;
		Element& neg(Element& r, const Element& a) const override;
		Element& inv(Element& r, const Element& a) const override;

		Element& mulin(Element& r, const Element& a) const override;
		Element& divin(Element& r, const Element& a) const override;
		Element& addin(Element& r, const Element& a) const override;
		Element& subin(Element& r, const Element& a) const override;
		Element& negin(Element& r) const override;
		Element& invin(Element& r) const override;

		// Functions defined in modular-mulprecomp
		//
		// void precomp_p (Compute_t& invp) const
		// Element& mul_precomp_p (Element& r, const Element& a, const Element& b, const Compute_t& invp) const
		//
		// void precomp_b (Compute_t& invb, const Element& b) const
		// void precomp_b (Compute_t& invb, const Element& b, const Compute_t& invp) const
		// Element& mul_precomp_b (Element& r, const Element& a, const Element& b, const Compute_t& invb) const

#include"modular-mulprecomp.inl"

		// -- axpy:   r <- a * x + y
		// -- axpyin: r <- a * x + r
		Element& axpy  (Element& r, const Element& a, const Element& x, const Element& y) const override;
		Element& axpyin(Element& r, const Element& a, const Element& x) const override;

		// -- axmy:   r <- a * x - y
		// -- axmyin: r <- a * x - r
		Element& axmy  (Element& r, const Element& a, const Element& x, const Element& y) const override;
		Element& axmyin(Element& r, const Element& a, const Element& x) const override;

		// -- maxpy:   r <- y - a * x
		// -- maxpyin: r <- r - a * x
		Element& maxpy  (Element& r, const Element& a, const Element& x, const Element& y) const override;
		Element& maxpyin(Element& r, const Element& a, const Element& x) const override;

		// ----- Random generators
		typedef ModularRandIter<Self_t> RandIter;
		typedef GeneralRingNonZeroRandIter<Self_t> NonZeroRandIter;
		template< class Random > Element& random(Random& g, Element& r) const
		{ return init(r, g()); }
		template< class Random > Element& nonzerorandom(Random& g, Element& a) const
		{ while (isZero(init(a, g())))
				;
			return a; }

		// --- IO methods
		std::istream& read (std::istream& s);
		std::ostream& write(std::ostream& s) const;
		std::istream& read (std::istream& s, Element& a) const;
		std::ostream& write(std::ostream& s, const Element& a) const;

	protected:
		// -- data representation of the domain:
		Residu_t _p;
		size_t _bitsizep;
	};
}

#include "givaro/modular-int32.inl"

#endif // __GIVARO_modular_int32_H
