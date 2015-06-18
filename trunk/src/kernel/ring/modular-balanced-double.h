// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
//          Clement Pernet <clement.pernet@gmail.com>
//          Brice Boyer <bboyer@imag.fr> (modified)
//          A. Breust (taken from FFLAS-FFPACK)
// ==========================================================================

#ifndef __GIVARO_modular_balanced_double_H
#define __GIVARO_modular_balanced_double_H

#include "givaro/givinteger.h"
#include "givaro/givcaster.h"
#include "givaro/givranditer.h"
#include "givaro/ring-interface.h"
#include "givaro/modular-general.h"

#include <iostream>

namespace Givaro
{
    template<class TAG> class ModularBalanced;

    template <>
    class ModularBalanced<double> : public virtual FiniteFieldInterface<double>
    {
    public:
	
	// ----- Exported types
	using Self_t = ModularBalanced<double>;
	enum { size_rep = sizeof(Element) };

	// ----- Constantes
	const Element zero = 0.0;
	const Element one  = 1.0;
	const Element mOne = -1.0;

	// ----- Constructors
	ModularBalanced()
	    : _p(0), _halfp(0), _mhalfp(0), _up(0.0)
	{}

	ModularBalanced(double p)
	    : _p(p)
	    , _halfp((_p - 1.f) / 2.f)
	    , _mhalfp(_halfp - _p + 1.f)
	    , _up(static_cast<uint32_t>(_p))
	{
	    assert(_p >= minCardinality());
	    assert(_p <= maxCardinality());
	}

	ModularBalanced(const Self_t& F)
	    : _p(F._p), _halfp(F._halfp), _mhalfp(F._mhalfp), _up(F._up)
	{}

	// ----- Accessors
	inline Element minElement() const override { return _mhalfp; }
	inline Element maxElement() const override { return _halfp; }

	// ----- Access to the modulus
	inline Element residu() const { return _p; }
	inline Element size() const { return _p; }
	inline Element characteristic() const { return _p; }
	inline Element cardinality() const { return _p; }
	template<class T> inline T& characteristic(T& p) const { return Caster(p,_p); }
	template<class T> inline T& cardinality(T& p) const { return Caster(p,_p); }
	
	static inline Element maxCardinality() { return 134217727; } // 2^12.5
	static inline Element minCardinality() { return 3.f; }

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
	    _halfp = F._halfp;
	    _mhalfp = F._mhalfp;
	    _up = F._up;
	    return *this;
	}

	// ----- Initialisation
	Element& init(Element& a) const;
	Element& init(Element& r, const float a) const;
	Element& init(Element& r, const double a) const;
	Element& init(Element& r, const int64_t a) const;
	Element& init(Element& r, const uint64_t a) const;
	Element& init(Element& r, const Integer& a) const;
	template<typename T> Element& init(Element& r, const T& a) const
	{ r = Caster<Element>(a); return reduce(r); }

	Element& assign(Element& r, const Element& a) const;

	// ----- Convert
	template<typename T> T& convert(T& r, const Element& a) const
	{ return r = static_cast<T>(a); }
	
	Element& reduce(Element& r, const Element& a) const;
	Element& reduce(Element& r) const;
	
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
	template< class Random > Element& random(const Random& g, Element& r) const
	{ return init(r, g()); }
	template< class Random > Element& nonzerorandom(const Random& g, Element& a) const
    	{ while (isZero(init(a, g())))
                ;
	    return a; }

	// --- IO methods
	std::ostream& write(std::ostream& s) const;
	std::istream& read (std::istream& s, Element& a) const;
	std::ostream& write(std::ostream& s, const Element& a) const;
	
    protected:

	Element _p;
	Element _halfp;
	Element _mhalfp;
	uint64_t _up;
    };
}

#include "givaro/modular-balanced-double.inl"

#endif 

