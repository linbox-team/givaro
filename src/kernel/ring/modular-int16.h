// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givzpz16std.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givzpz16std.h,v 1.16 2011-02-02 17:16:43 bboyer Exp $
// ==========================================================================
//
//  Modified by Pascal Giorgi on 2002/02/13  (pascal.giorgi@ens-lyon.fr)
//  Modified by Alexis Breust on 2015/01/06  (alexis.breust@imag.fr)

/*! @file givzpz16std.h
 * @ingroup zpz
 * @brief   Arithmetic on Z/pZ, with p a prime number less than 2^14.
 *   Modulo typedef is a signed long number. In case it was modified
 *   then Bézout algorithm must be changed (coefficient can be negative).
 */
 
#ifndef __GIVARO_zpz16std_H
#define __GIVARO_zpz16std_H

#include "givaro/givbasictype.h"
#include "givaro/giverror.h"
#include "givaro/givranditer.h"
#include "givaro/ring-interface.h"

namespace Givaro {

template<class TAG> class Modular;

/*! @brief This class implement the standard arithmetic with Modulo Elements.
 * - The representation of an integer a in Zpz is the value a % p
 * - m max is 32768
 * - p max is 32749
 * .
 */
template<>
class Modular<int16_t> : public RingInterface<int16_t>
{
public:

	// ----- Exported Types and constantes
	typedef Modular<int16_t> Self_t;
	typedef uint16_t Residu_t;
	typedef uint32_t Compute_t;
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
	inline Element minElement() const override { return zero; }
	inline Element maxElement() const override { return mOne; }

	// ----- Access to the modulus
	inline Residu_t residu() const { return _p; }
	inline Residu_t size() const { return _p; }
	inline Residu_t characteristic() const { return _p; }
	inline Residu_t cardinality() const { return _p; }
	template<class T> inline T& characteristic(T& p) const { return p = _p; }
	template<class T> inline T& cardinality(T& p) const { return p = _p; }
	static inline Residu_t getMaxModulus() { return 32768; } // 2^15 (because all computations are done modulo 2^32)
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
		return *this;
	}
	
	// ----- Initialisation
	Element& init(Element& a) const;
	Element& init(Element& r, const long a) const ;
	Element& init(Element& r, const unsigned long a) const ;
	Element& init(Element& a, const int i) const ;
	Element& init(Element& a, const unsigned int i) const ;
	Element& init(Element& a, const Integer& i) const ;
	Element& init(Element& a, const double i) const;
	Element& init(Element& a, const float i) const;
	void init(const size_t, Array a, constArray b) const;
	
	Element& assign(Element& r, const Element a) const;
	void assign(const size_t sz, Array r, constArray a ) const;

	// ----- Convert
	float& convert(float& r, const Element a ) const { return r = (float)a ;}
	double& convert(double& r, const Element a ) const { return r = (double)a ;}
	long int& convert(long int& r, const Element a) const { return r = (long int)a;}
	unsigned long int& convert(unsigned long int& r, const Element a) const { return r = (unsigned long int)a;}
	int16_t& convert(int16_t& r, const Element a ) const { return r = (int16_t)a ;}
	int32_t& convert(int32_t& r, const Element a ) const { return r = (int32_t)a ;}
	uint16_t& convert(uint16_t& r, const Element a ) const { return r = (uint16_t)a ;}
	uint32_t& convert(uint32_t& r, const Element a ) const { return r = (uint32_t)a ;}
	Integer& convert(Integer& i, const Element a) const { unsigned long ur; return i = (Integer)convert(ur, a);	}

	inline Element& reduce (Element& x, const Element& y) const { return x = Element(y % _p); }
	inline Element& reduce (Element& x) const { return x %= _p; }
	
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
    template< class Random > Element& random(const Random& g, Element& r) const	{ return init(r, g()); }
    template< class Random > Element& nonzerorandom(const Random& g, Element& a) const
    	{ while (isZero(init(a, g())));
    	  return a; }
	
	// --- IO methods
	std::istream& read (std::istream& s);
	std::ostream& write(std::ostream& s) const;
	std::istream& read (std::istream& s, Element& a) const;
	std::ostream& write(std::ostream& s, const Element a) const;
	
protected:
	// -- Modular inverse, d = a*u + b*v
	int32_t& gcdext(int32_t& d, int32_t& u, int32_t& v, const int32_t a, const int32_t b) const;
	int32_t& invext(int32_t& u, const int32_t a, const int32_t b) const;

protected:
	// -- data representation of the domain:
	Residu_t _p;
};

} // namespace Givaro

#include "givaro/modular-int16.inl"

#endif // __GIVARO_zpz16std_H
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
