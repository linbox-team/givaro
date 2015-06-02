// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: A. Breust <alexis.breust@imag.fr>
// ==========================================================================

#ifndef __GIVARO_montgomery_ruint_H
#define __GIVARO_montgomery_ruint_H

#include "recint/ruint.h"
#include "recint/rmgmodule.h"
#include "givaro/ring-interface.h"

namespace Givaro
{

//! @brief The recint-based Montgomery ring.
//! You can only odd moduli.
//! An integer (a mod p) is stored as (a * r mod 2^{2^K}) with (r = 2^{2^K} mod p).

template<size_t K>
class Montgomery<RecInt::ruint<K>> : public RingInterface<RecInt::ruint<K>>
{
public:

	// ----- Exported Types and constantes
	using Element = RecInt::ruint<K>;
	using LargeElement = RecInt::ruint<K+1>;
	using Self_t = Montgomery<RecInt::ruint<K>>;
	using Residu_t = RecInt::ruint<K>;
	enum { size_rep = sizeof(Residu_t) };

	// ----- Representation of vector of the Element
	typedef Element* Array;

	// ----- Constantes
	const Element zero;
	const Element one;
	const Element mOne;

	// ----- Constructors
	Montgomery()
	    :  zero(0), one(1), mOne(-1)
        , _p(0), _p1(0), _r(0), _r3(0)
    {}

	Montgomery(const Residu_t& p)
	    : zero(0)
	    , _p(p)
	{
        RecInt::arazi_qi(_p1, -_p); // p1 = -inv(p) mod 2^(2^K)
        RecInt::mod_n(_r, -_p, _p); // r = 2^(2^K) mod p

        RecInt::mul(_r3, _r, _r);
        RecInt::mul(_r3, _r);       // r3 = r^3 mod p

        RecInt::copy(const_cast<Element&>(one), _r);
        to_mg(const_cast<Element&>(mOne), _p - 1u);

        assert(_p & 1u != 0u);
	    assert(_p >= getMinModulus());
	    assert(_p <= getMaxModulus());
	}

	Montgomery(const Self_t& F)
	    : zero(F.zero), one(F.one), mOne(F.mOne)
	    , _p(F._p), _p1(F._p1), _r(F._r)
	{}

	// ----- Accessors
	inline Element minElement() const final { return zero; }
	inline Element maxElement() const final { return mOne; }

	// ----- Access to the modulus
	inline Residu_t residu() const { return _p; }
	inline Residu_t size() const { return _p; }
	inline Residu_t characteristic() const { return _p; }
	inline Residu_t cardinality() const { return _p; }
	template<class T> inline T& characteristic(T& p) const { return p = _p; }
	template<class T> inline T& cardinality(T& p) const { return p = _p; }

	static inline Residu_t getMaxModulus() { return -1; }
	static inline Residu_t getMinModulus() { return 2; }

	// ----- Checkers
	inline bool isZero(const Element& a) const final { return a == zero; }
	inline bool isOne (const Element& a) const final { return a == one; }
	inline bool isMOne(const Element& a) const final { return a == mOne; }
	inline bool areEqual(const Element& a, const Element& b) const final { return a == b; }
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
		return *this;
	}

	// ----- Initialisation
	Element& init(Element& r) const;
	Element& init(Element& r, const double a)   const;
	Element& init(Element& r, const Integer& a) const;

	Element& assign(Element& r, const Element& a) const;

	// ----- Convert
	Integer& convert(Integer& i, const Element a) const { unsigned long ur; return i = (Integer)convert(ur, a);	}
	template<typename XXX> XXX& convert(XXX& r, const Element a) const { return r = static_cast<XXX>(a) ;}

	inline Element& reduce (Element& x, const Element& y) const { return x = Element(y % _p); }
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
    template< class Random > Element& random(const Random& g, Element& r) const	{ return init(r, g()); }
    template< class Random > Element& nonzerorandom(const Random& g, Element& a) const
    	{ while (isZero(init(a, g())));
    	  return a; }

	// --- IO methods
	std::istream& read (std::istream& s);
	std::ostream& write(std::ostream& s) const;
	std::istream& read (std::istream& s, Element& a) const;
	std::ostream& write(std::ostream& s, const Element& a) const;

protected:

    // Internal montgomery-reduction. a <- b * r^{-1}
    Element& mg_reduc(Element& a, const Element& b) const;
    Element& mg_reduc(Element& a, const LargeElement& b) const;

    // a <- b * r
    Element& to_mg(Element& a, const Element& b) const;
    Element& to_mg(Element& a) const;

protected:

    // p is the module (must be odd and > 1)
    RecInt::ruint<K> _p;
    // p1 = -inv(p) mod 2^(2^K)
    RecInt::ruint<K> _p1;
    // r = 2^(2^K) mod p
    RecInt::ruint<K> _r;
    // r3 = r^3 mod p - used to compute inverse
    RecInt::ruint<K> _r3;
};

} // namespace Givaro

#include "givaro/montgomery-ruint.inl"

#endif
