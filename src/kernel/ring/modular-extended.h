/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
// ==========================================================================
// Copyright(c)'1994-2016 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: Bastien Vialla <bastien.vialla@lirmm.fr>
// ==========================================================================

#ifndef __GIVARO_MODULAR_EXTENDED_H
#define __GIVARO_MODULAR_EXTENDED_H

#include "givaro/givconfig.h"

#include "givaro/givranditer.h"
#include "givaro/ring-interface.h"
#include "givaro/modular-general.h"

namespace Givaro{
	/*
	 *
	 * Modular double/float allowing big moduli
	 * !!: RandIter does not works, use your own random
	 *
	 */
	template<class _Element>
	class ModularExtended : public virtual FiniteFieldInterface<_Element>
	{
	public:

		typedef _Element Element;
		typedef Element* Element_ptr ;
		typedef const Element ConstElement;
		typedef const Element* ConstElement_ptr;
		// ----- Exported Types and constantes
		typedef ModularExtended<_Element> Self_t;
		using Compute_t = _Element;
		typedef uint64_t Residu_t;
		enum { size_rep = sizeof(Residu_t) };

	private:
		// Verkampt Split
		inline void split(const Element x, Element &x_h, Element &x_l) const {
			Element c;
			if(std::is_same<Element, double>::value){
				c = (Element)((1 << 27)+1);
			}else if(std::is_same<Element, float>::value){
				c = (Element)((1 << 13)+1);
			}
			c = c*x;
			x_h = c-(c-x);
			x_l = x - x_h;
		}

		// Dekker mult, a * b = s + t
		inline void mult_dekker(const Element a, const Element b, Element &s, Element &t) const{
			s = a*b;
			Element ah, al, bh, bl;
			split(a, ah, al);
			split(b, bh, bl);
			t = (al*bl-(((s-ah*bh)-al*bh)-ah*bl));
		}

	public:
		// ----- Constantes
		const Element zero = 0.0;
		const Element one = 1.0;
		const Element mOne = -1.0;

		// ----- Constructors
		ModularExtended() = default;

		template<class XXX> ModularExtended(const XXX& p)
			: zero(0.0), one(1.0), mOne((Element)p - 1.0), _p((Element)p), _invp((Element)1/(Element)_p), _negp(-_p), _lp((Residu_t)p)
		{
			assert(_p >= minCardinality());
			assert(_p <= maxCardinality());
		}

		// ----- Accessors
		inline Element minElement() const override { return zero; }
		inline Element maxElement() const override { return mOne; }

		// ----- Access to the modulus
		inline Residu_t residu() const { return _lp; }
		inline Residu_t size() const { return _lp; }
		inline Residu_t characteristic() const { return _lp; }
		template<class T> inline T& characteristic(T& p) const { return p = _lp; }
		inline Residu_t cardinality() const { return _lp; }
		template<class T> inline T& cardinality(T& p) const { return p = _lp; }
		static inline Residu_t maxCardinality() {
			if(std::is_same<Element, double>::value)
				return 1125899906842623; // 2^(52-2) - 1
			else if(std::is_same<Element, float>::value)
				return 2097151; // 2^(23-2) - 1
		}
		static inline Residu_t minCardinality() { return 2; }

		// ----- Checkers
		inline bool isZero(const Element& a) const override { return a == zero; }
		inline bool isOne (const Element& a) const override { return a == one; }
		inline bool isMOne(const Element& a) const override { return a == mOne; }
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
			_negp = F._negp;
			_invp = F._invp;
			_lp= F._lp;
			return *this;
		}



		// ----- Initialisation
		inline Element& init (Element& x) const
		{
			return x = zero;
		}

		template<typename T>
		Element& init (Element& r, T a) const{
			r = Caster<Element>(a);
			return reduce(r);
		}

		Element &assign (Element &x, const Element &y) const{
			return x = y;
		}

		// ----- Convert and reduce
		template<typename T> T& convert(T& r, const Element& a) const
		{ return r = static_cast<T>(a); }

		Element& reduce (Element& x, const Element& y) const{
			x=y;
			return reduce(x);
		}

		Element& reduce (Element& x) const ;

		// ----- Classic arithmetic
		Element& mul(Element& r, const Element& a, const Element& b) const override;

		Element& div(Element& r, const Element& a, const Element& b) const override{
			return mulin(inv(r, a), b);
		}
		Element& add(Element& r, const Element& a, const Element& b) const override {
			r = a + b;
			if(r >= _p)
				r += _negp;
			return r;
		}
		Element& sub(Element& r, const Element& a, const Element& b) const override {
			r = a - b;
			if(r < 0)
				r += _p;
			return r;
		}
		Element& neg(Element& r, const Element& a) const override {
			r = -a;
			if(r < 0)
				r += _p;
			return r;
		}
		Element& inv(Element& x, const Element& y) const override{
			if(std::is_same<Element, double>::value){
				int64_t x_int, y_int, tx, ty;
				x_int = int64_t(_lp);
				y_int = int64_t(y);
				tx = 0;
				ty = 1;

				while (y_int != 0) {
					// always: gcd (modulus,residue) = gcd (x_int,y_int)
					//         sx*modulus + tx*residue = x_int
					//         sy*modulus + ty*residue = y_int
					int64_t q = x_int / y_int; // integer quotient
					int64_t temp = y_int;  y_int  = x_int  - q * y_int;
					x_int  = temp;
					temp = ty; ty = tx - q * ty;
					tx = temp;
				}

				if (tx < 0) tx += int64_t(_p);

				// now x_int = gcd (modulus,residue)
				return x = Element(tx);
			}else if(std::is_same<Element, float>::value){
				int32_t x_int, y_int, tx, ty;
				x_int = int32_t(_lp);
				y_int = int32_t(y);
				tx = 0;
				ty = 1;

				while (y_int != 0) {
					// always: gcd (modulus,residue) = gcd (x_int,y_int)
					//         sx*modulus + tx*residue = x_int
					//         sy*modulus + ty*residue = y_int
					int32_t q = x_int / y_int; // integer quotient
					int32_t temp = y_int;  y_int  = x_int  - q * y_int;
					x_int  = temp;
					temp = ty; ty = tx - q * ty;
					tx = temp;
				}

				if (tx < 0) tx += int32_t(_p);

				// now x_int = gcd (modulus,residue)
				return x = Element(tx);
			}
		}

		Element& mulin(Element& r, const Element& a) const override {
			return mul(r, r, a);
		}
		Element& divin(Element& r, const Element& y) const override{
			Element iy;
			return mulin(r, inv(iy, y));
		}
		Element& addin(Element& r, const Element& a) const override {
			return add(r, r, a);
		}
		Element& subin(Element& r, const Element& a) const override {
			return sub(r, r, a);
		}
		Element& negin(Element& r) const override {
			return neg(r, r);
		}
		Element& invin(Element& r) const override {
			return inv(r, r);
		}

		// -- axpy:   r <- a * x + y
		// -- axpyin: r <- a * x + r
		Element& axpy  (Element& r, const Element& a, const Element& x, const Element& y) const override {
			Element tmp;
			mul(tmp, a, x);
			return add(r, tmp, y);
		}
		Element& axpyin(Element& r, const Element& a, const Element& x) const override {
			Element tmp(r);
			return axpy(r, a, x, tmp);
		}

		// -- axmy:   r <- a * x - y
		// -- axmyin: r <- a * x - r
		Element& axmy  (Element& r, const Element& a, const Element& x, const Element& y) const override {
			Element tmp;
			mul(tmp, a, x);
			return sub(r, tmp, y);
		}
		Element& axmyin(Element& r, const Element& a, const Element& x) const override {
			return axmy(r, a, x, r);
		}

		// -- maxpy:   r <- y - a * x
		// -- maxpyin: r <- r - a * x
		Element& maxpy  (Element& r, const Element& a, const Element& x, const Element& y) const override {
			Element tmp;
			mul(tmp, a, x);
			return sub(r, y, tmp);
		}
		Element& maxpyin(Element& r, const Element& a, const Element& x) const override {
			return maxpy(r, a, x, r);
		}

		// ----- Random generators
		typedef ModularRandIter<Self_t> RandIter;
		typedef GeneralRingNonZeroRandIter<Self_t> NonZeroRandIter;
		template< class Random > Element& random(const Random& g, Element& r) const { return init(r, g()); }
		template< class Random > Element& nonzerorandom(const Random& g, Element& a) const {
			while (isZero(init(a, g())));
			return a;
		}

		// --- IO methods
		std::istream& read (std::istream& s);
		std::ostream& write(std::ostream& s) const;
		std::istream& read (std::istream& s, Element& a) const;
		std::ostream& write(std::ostream& s, const Element& a) const;

	protected:
		Element _p = 0;
		Element _invp = 0;
		Element _negp = 0;
		Residu_t _lp = 0;

	};

	// ----------------
	// ----- IO methods

	template<>
	inline
	std::ostream &ModularExtended<float>::write (std::ostream &os) const
	{
		return os << "ModularExtended<float> mod " << _p;
	}

	template<>
	inline
	std::ostream &ModularExtended<double>::write (std::ostream &os) const
	{
		return os << "ModularExtended<double> mod " << _p;
	}

	template<typename _Element>
	inline
	std::istream &ModularExtended<_Element>::read (std::istream &is)
	{
		is >> _p;
		return is;
	}

	template<typename _Element>
	inline
	std::ostream &ModularExtended<_Element>::write (std::ostream &os, const Element &x) const
	{
		return os << static_cast<int64_t>(x);
	}

	template<typename _Element>
	inline
	std::istream &ModularExtended<_Element>::read (std::istream &is, Element &x) const
	{
		int64_t tmp;
		is >> tmp;
		init(x,tmp);
		return is;
	}

	//Element &init (Element &x) const{
	//	return x = zero;
	//}

	//template<class XXX> Element& init(Element & x, const XXX & y) const{
	//	// TODO : Bad init, take the code from modular<double>
	//	x=Element(y);
	//	return reduce(x);
	//}

	// --------------------
	// ----- Initialisation de Modular<double>

	template<>
	template<>
	typename ModularExtended<double>::Element&
	ModularExtended<double>::init (typename ModularExtended<double>::Element& x, const int64_t y) const
	{
		x = static_cast<Element>(std::abs(y) % _lp);
		if (y < 0) negin(x);
		return x;
	}

	template<>
	template<>
	typename ModularExtended<double>::Element&
	ModularExtended<double>::init (typename ModularExtended<double>::Element& x, const uint64_t y) const
	{
		return x = static_cast<Element>(y % (uint64_t)(_lp));
	}

	template<>
	template<>
	typename ModularExtended<double>::Element&
	ModularExtended<double>::init (typename ModularExtended<double>::Element& x, const Integer& y) const
	{
		x = static_cast<Element>(y % _lp);
		if (x < 0) x += _p;
		return x;
	}

	template<>
	template<>
	typename ModularExtended<double>::Element&
	ModularExtended<double>::init (typename ModularExtended<double>::Element& x, const Element& y) const
	{
		return x = y;
	}

	// --------------------
	// ----- Initialisation de ModularExtended<float>

	template<>
	template<>
	typename ModularExtended<float>::Element&
	ModularExtended<float>::init(typename ModularExtended<float>::Element& r, const double a) const
	{
		r = static_cast<float>(std::fmod(a, _p));
		if (r < 0.f) r += _p;
		return r;
	}

	template<>
	template<>
	typename ModularExtended<float>::Element&
	ModularExtended<float>::init(typename ModularExtended<float>::Element& r, const int32_t a) const
	{
		r = static_cast<Element>(std::abs(a) % _lp);
		if (a < 0) negin(r);
		return r;
	}

	template<>
	template<>
	typename ModularExtended<float>::Element&
	ModularExtended<float>::init(typename ModularExtended<float>::Element& r, const uint32_t a) const
	{
		return r = static_cast<Element>(a % uint32_t(_lp));
	}

	template<>
	template<>
	typename ModularExtended<float>::Element&
	ModularExtended<float>::init(typename ModularExtended<float>::Element& r, const int64_t a) const
	{
		r = static_cast<Element>(std::abs(a) % int64_t(_lp));
		if (a < 0) negin(r);
		return r;
	}

	template<>
	template<>
	typename ModularExtended<float>::Element&
	ModularExtended<float>::init(typename ModularExtended<float>::Element& r, const uint64_t a) const
	{
		return r = static_cast<Element>(a % uint64_t(_lp));
	}

	template<>
	template<>
	typename ModularExtended<float>::Element&
	ModularExtended<float>::init(typename ModularExtended<float>::Element& r, const Integer& a) const
	{
		r = static_cast<Element>(a % _lp);
		if (a < 0) negin(r);
		return r;
	}

	// --------------
	// Multiplication
	// --------------
	template<>
	typename ModularExtended<float>::Element&
	ModularExtended<float>::mul(typename ModularExtended<float>::Element& r,
								const typename ModularExtended<float>::Element& a,
								const typename ModularExtended<float>::Element& b) const {
		Element abh, abl, pql, pqh, q;
#ifdef FP_FAST_FMAF
		abh = a * b;
		abl = fma(a, b, -abh);
		q = std::floor(abh*_invp);
		pql = fma (-q, _p, abh);
		r = abl + pql;
#else
		mult_dekker(a, b, abh, abl);
		q = std::floor(abh*_invp);
		mult_dekker(-q, _p, pqh, pql);
		r = (abh + pqh) + (abl + pql);
#endif
		if(r >= _p)
			r-= _p;
		else if(r < 0)
			r += _p;
#ifndef NDEBUG
		assert((r < _p) && (r >= 0));
#endif
		return r;
	}

	template<>
	typename ModularExtended<double>::Element&
	ModularExtended<double>::mul(typename ModularExtended<double>::Element& r,
								const typename ModularExtended<double>::Element& a,
								const typename ModularExtended<double>::Element& b) const {
		Element abh, abl, pql, pqh, q;
#ifdef FP_FAST_FMA
		abh = a * b;
		abl = fma(a, b, -abh);
		q = std::floor(abh*_invp);
		pql = fma (-q, _p, abh);
		r = abl + pql;
#else
		mult_dekker(a, b, abh, abl);
		q = std::floor(abh*_invp);
		mult_dekker(-q, _p, pqh, pql);
		r = (abh + pqh) + (abl + pql);
#endif
		if(r >= _p)
			r-= _p;
		else if(r < 0)
			r += _p;
#ifndef NDEBUG
		assert((r < _p) && (r >= 0));
#endif
		return r;
	}

	// --------------
	//    Reduction
	// --------------
	template<>
	typename ModularExtended<float>::Element&
	ModularExtended<float>::reduce (typename ModularExtended<float>::Element& a) const{
		Element pql, pqh, q;
#ifdef FP_FAST_FMAF
		q = std::floor(a*_invp);
		pql = fma (-q, _p, a);
		a = pql;
#else
		q = std::floor(a*_invp);
		mult_dekker(-q, _p, pqh, pql);
		a = (a + pqh) + pql;
#endif
		if(a >= _p)
			a-= _p;
		else if(a < 0)
			a += _p;
#ifndef NDEBUG
		assert((a < _p) && (a >= 0));
#endif
		return a;
	}

	template<>
	typename ModularExtended<double>::Element&
	ModularExtended<double>::reduce (typename ModularExtended<double>::Element& a) const{
		Element pql, pqh, q;
#ifdef FP_FAST_FMA
		q = std::floor(a*_invp);
		pql = fma (-q, _p, a);
		a = pql;
#else
		q = std::floor(a*_invp);
		mult_dekker(-q, _p, pqh, pql);
		a = (a + pqh) + pql;
#endif
		if(a >= _p)
			a-= _p;
		else if(a < 0)
			a += _p;
#ifndef NDEBUG
		assert((a < _p) && (a >= 0));
#endif
		return a;
	}


} // Givaro

#endif //__GIVARO_MODULAR_EXTENDED_H
