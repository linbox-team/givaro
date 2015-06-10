// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: Clement Pernet <clement.pernet@gmail.com>
//          Pascal Giorgi  <pascal.giorgi@ens-lyon.fr>
//          A. Breust (taken from FFLAS-FFPACK)
// ==========================================================================

#ifndef __GIVARO_modular_double_H
#define __GIVARO_modular_double_H

#include <float.h>

#include "givaro/givranditer.h"
#include "givaro/ring-interface.h"
#include "givaro/modular-general.h"

namespace Givaro
{

template <>
class Modular<double, double> : public RingInterface<double>
{
public:
	// ----- Exported Types and constantes
	typedef Modular<double> Self_t;
	typedef uint64_t Residu_t;
	enum { size_rep = sizeof(Residu_t) };

	// ----- Constantes
	const Element zero;
	const Element one;
	const Element mOne;

	// ----- Constructors
	Modular()
	: zero(0.0), one(1.0), mOne(-1.0), _p(0.0), _lp(0) {}

	template<class XXX> Modular(const XXX& p)
	: zero(0.0), one(1.0), mOne((Element)p - 1.0), _p((double)p), _lp((Residu_t)p)
	{
	    assert(_p >= getMinModulus());
	    assert(_p <= getMaxModulus());
	}

	Modular(const Self_t& F)
	: zero(F.zero), one(F.one), mOne(F.mOne), _p(F._p), _lp(F._lp) {}

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
	static inline Residu_t getMaxModulus() { return 67108864; }
	static inline Residu_t getMinModulus() { return 2; }

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
		_lp= F._lp;
		return *this;
	}

	// ----- Initialisation
	Element& init (Element& x) const;
	Element& init (Element& x, const int64_t y) const;
	Element& init (Element& x, const uint64_t y) const;
	Element& init (Element& x, const Integer& y) const;
	template<typename T> Element& init(Element& r, const T& a) const
	{ r = Caster<Element>(a); return reduce(r); }

	Element& assign (Element& x, const Element& y) const;

	// ----- Convert and reduce
	template<typename T> T& convert(T& r, const Element& a) const
	{ return r = static_cast<T>(a); }

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
	double _p;
	Residu_t _lp;

};

} // Givaro

#include "givaro/modular-double.inl"

#endif

