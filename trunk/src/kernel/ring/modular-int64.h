// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givzpz64std.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givzpz64std.h,v 1.21 2011-02-04 14:11:46 jgdumas Exp $
// ==========================================================================

/*! @file givzpz64std.h
 * @ingroup zpz
 * @brief Zpz on 64bit words
 *   Arithmetic on Z/pZ, with p a prime number less than 2^64
 *   Modulo typedef is a signed long number. In case it was modified
 *   then Bézout algorithm must be changed (coefficient can be negative).
 */

#ifndef __GIVARO_zpz64std_H
#define __GIVARO_zpz64std_H

#include "givaro/givinteger.h"
#include "givaro/givcaster.h"
#include "givaro/givranditer.h"
#include "givaro/ring-interface.h"
#include "givaro/modular-general.h"

namespace Givaro {

/*! @brief This class implement the standard arithmetic with Modulo Elements.
 * - The representation of an integer a in Zpz is the value a % p
 * - m max is 3037000499
 * - p max is 3037000493
 * .
 */
template<typename COMP>
class Modular<int64_t, COMP> : public RingInterface<int64_t>
{
public:

	// ----- Exported Types and constantes
	using Self_t = Modular<int64_t, COMP>;
	using Residu_t = uint64_t;
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
	: zero(0), one(1), mOne(-1), _p(0) {}

	Modular(Residu_t p)
	: zero(0), one(1), mOne((Element)p-1), _p(p)
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
	inline Residu_t cardinality() const { return _p; }
	template<class T> inline T& characteristic(T& p) const { return p = _p; }
	template<class T> inline T& cardinality(T& p) const { return p = _p; }
	static inline Residu_t getMaxModulus() { return 3037000499ULL; } 
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
	Element& init (Element& x) const;
	Element& init (Element& x, const Integer& y) const;
	template<typename T> Element& init(Element& r, const T& a) const
	{ r = Caster<Element>(a); return reduce(r); }
	void init(const size_t, Array a, constArray b) const;

	Element& assign (Element& x, const Element& y) const;
	void assign(const size_t sz, Array r, constArray a ) const;
    
	// ----- Convert and reduce
	template<typename T> T& convert(T& r, const Element& a) const
	{ return r = static_cast<T>(a); }

	Element& reduce (Element& x, const Element& y) const;
	Element& reduce (Element& x) const;
	
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
	
	// ----- Classic arithmetic on arrays
	void mul(const size_t sz, Array r, constArray a, constArray b) const;
	void mul(const size_t sz, Array r, constArray a, Element b) const;
	void div(const size_t sz, Array r, constArray a, constArray b) const;
	void div(const size_t sz, Array r, constArray a, Element b) const;
	void add(const size_t sz, Array r, constArray a, constArray b) const;
	void add(const size_t sz, Array r, constArray a, Element b) const;
	void sub(const size_t sz, Array r, constArray a, constArray b) const;
	void sub(const size_t sz, Array r, constArray a, Element b) const;
	void neg(const size_t sz, Array r, constArray a) const;
	void inv(const size_t sz, Array r, constArray a) const;

	void axpy (const size_t sz, Array r, constArray a, constArray x, constArray c) const;
	void axpyin (const size_t sz, Array r, constArray a, constArray x) const;
	void axmy (const size_t sz, Array r, constArray a, constArray x, constArray c) const;
	void maxpyin (const size_t sz, Array r, constArray a, constArray x) const;
	
	// <- \sum_i a[i], return 1 if a.size() ==0,
	Element& reduceadd ( Element& r, const size_t sz, constArray a ) const;

	// <- \prod_i a[i], return 1 if a.size() ==0,
	Element& reducemul ( Element& r, const size_t sz, constArray a ) const;

	// <- \sum_i a[i] * b[i]
	Element& dotprod ( Element& r, const size_t sz, constArray a, constArray b ) const;
	Element& dotprod ( Element& r, const int bound, const size_t sz, constArray a, constArray b ) const;

	// a -> r: uint32_t to double
	void i2d ( const size_t sz, double* r, constArray a ) const;

	// a -> r % p: double to uint32_t % p
	void d2i ( const size_t sz, Array r, const double* a ) const;

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
	std::ostream& write(std::ostream& s, const Element a) const;

protected:
	// -- data representation of the domain:
	Residu_t _p;
};

} // namespace Givaro

#include "givaro/modular-int64.inl"
#endif // __GIVARO_zpz64std_H

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
