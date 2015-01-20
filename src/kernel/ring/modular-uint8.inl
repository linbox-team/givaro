// ==========================================================================
// Copyright(c)'1994-2014 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: A. Breust (taken from FFLAS-FFPACK)
// ==========================================================================

#ifndef __GIVARO_modular_uint8_INL
#define __GIVARO_modular_uint8_INL

#include <cmath> // fmod

#include "givaro/modular-defines.h"

namespace Givaro {

	// --------------------
	// ----- Initialisation
	
	inline Modular<uint8_t>::Element &Modular<uint8_t>::init (Element &x) const
	{
		return x = zero ;
	}

	inline Modular<uint8_t>::Element &Modular<uint8_t>::init (Element &x, const int32_t &y) const
	{
		x = (Element)(std::abs(y) % (uint64_t)(_p));
		if (y < 0) x = _p - x;
		return x;
	}

	inline Modular<uint8_t>::Element &Modular<uint8_t>::init (Element &x, const int64_t &y) const
	{
		x = (Element)(std::abs(y) % (uint64_t)(_p));
		if (y < 0) x = _p - x;
		return x;
	}

	inline Modular<uint8_t>::Element &Modular<uint8_t>::init(Element &x, const uint32_t &y) const
	{
		return x = (Element)( y >= (uint64_t)_p ? y % (uint64_t)(_p) : y);
	}

	inline Modular<uint8_t>::Element &Modular<uint8_t>::init (Element &x, const uint64_t &y) const
	{
		return x = (Element)( y >= (uint64_t)_p ? y % (uint64_t)(_p) : y);
	}

	inline Modular<uint8_t>::Element &Modular<uint8_t>::init (Element &x, const double &y) const
	{
		double z = fmod(y, (double)_p);
		if (z < 0) z += (double) _p;
		return x = (Element) (z);
	}

	template <class XXX>
	inline Modular<uint8_t>::Element &Modular<uint8_t>::init(Element &x, const XXX &y) const
	{
		return init(x, double(y));
	}
	
	inline Modular<uint8_t>::Element &Modular<uint8_t>::assign (Element &x, const Element &y) const
	{
		return x = y;
	}
	
	// -------------
	// ----- Convert

	inline double &Modular<uint8_t>::convert (double &x, const Element &y) const
	{
		return x = (double)y;
	}

	inline float &Modular<uint8_t>::convert (float &x, const Element &y) const
	{
		return x = (float)y;
	}

	// ------------------------
	// ----- Classic arithmetic
	
	inline Modular<uint8_t>::Element &Modular<uint8_t>::add
		(Element &x, const Element &y, const Element &z) const
	{
		__GIVARO_MODULAR_INTEGER_ADD(x,_p,y,z);
		return x;
	}

	inline Modular<uint8_t>::Element &Modular<uint8_t>::sub
		(Element &x, const Element &y, const Element &z) const
	{
		return __GIVARO_MODULAR_INTEGER_SUB(x,_p,y,z);
	}

	inline Modular<uint8_t>::Element &Modular<uint8_t>::mul
		(Element &x, const Element &y, const Element &z) const
	{
		return __GIVARO_MODULAR_INTEGER_MUL(x,_p,y,z);
	}

	inline Modular<uint8_t>::Element &Modular<uint8_t>::div
		(Element &x, const Element &y, const Element &z) const
	{
		return mulin(inv(x, z), y);
	}

	inline Modular<uint8_t>::Element &Modular<uint8_t>::neg
		(Element &x, const Element &y) const
	{
		return __GIVARO_MODULAR_INTEGER_NEG(x,_p,y);
	}

	inline Modular<uint8_t>::Element &Modular<uint8_t>::inv
		(Element &x, const Element &y) const
	{
		// The extended Euclidean algorithm
		int64_t x_int, y_int, tx, ty;
		x_int = _p;
		y_int = y;
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

		if (tx < 0) tx += _p;

		// now x_int = gcd (modulus,residue)
		return x = Element(tx);
	}

	inline Modular<uint8_t>::Element &Modular<uint8_t>::addin
		(Element &x, const Element &y) const
	{
		__GIVARO_MODULAR_INTEGER_ADDIN(x,_p,y);
		return x;
	}

	inline Modular<uint8_t>::Element &Modular<uint8_t>::subin
		(Element &x, const Element &y) const
	{
		__GIVARO_MODULAR_INTEGER_SUBIN(x,_p,y);
		return x;
	}

	inline Modular<uint8_t>::Element &Modular<uint8_t>::mulin
		(Element &x, const Element &y) const
	{
		return __GIVARO_MODULAR_INTEGER_MULIN(x,_p,y);
	}

	inline Modular<uint8_t>::Element &Modular<uint8_t>::divin
		(Element &x, const Element &y) const
	{
		Modular<uint8_t>::Element iy;
		return mulin(x, inv(iy, y));
	}

	inline Modular<uint8_t>::Element &Modular<uint8_t>::negin
		(Element &x) const
	{
		return __GIVARO_MODULAR_INTEGER_NEGIN(x,_p);
	}

	inline Modular<uint8_t>::Element &Modular<uint8_t>::invin
		(Element &x) const
	{
		return inv(x, x);
	}

	// -- axpy: r <- a * x + y
	inline Modular<uint8_t>::Element &Modular<uint8_t>::axpy
		(Element &r, const Element &a, const Element &x, const Element &y) const
	{
		__GIVARO_MODULAR_INTEGER_MULADD(r, _p, a, x, y);
		return r;
	}

	inline Modular<uint8_t>::Element &Modular<uint8_t>::axpyin
		(Element &r, const Element &a, const Element &x) const
	{
		__GIVARO_MODULAR_INTEGER_MULADDIN(r, _p, a, x);
		return r;
	}
	
	// -- axmy: r <- a * x - y
	inline Modular<uint8_t>::Element &Modular<uint8_t>::axmy
		(Element& r, const Element &a, const Element &x, const Element &y) const
	{
		__GIVARO_MODULAR_INTEGER_MULSUB(r, _p, a, x, y);
		return r;
	}
	
	inline Modular<uint8_t>::Element &Modular<uint8_t>::axmyin
		(Element& r, const Element &a, const Element &x) const
	{
		maxpyin(r,a,x);
		return negin(r);
	}
	
	// -- maxpy:   r <- y - a * x
	inline Modular<uint8_t>::Element& Modular<uint8_t>::maxpy
		(Element& r, const Element& a, const Element& x, const Element& y) const
	{
		r = y;
		__GIVARO_MODULAR_INTEGER_SUBMULIN(r, _p, a, x);
		return r;
	}
		
	inline Modular<uint8_t>::Element& Modular<uint8_t>::maxpyin
		(Element& r, const Element& a, const Element& x) const
	{
		__GIVARO_MODULAR_INTEGER_SUBMULIN(r, _p, a, x);
		return r;
	}
	
	// ----------------
	// ----- IO methods

	std::ostream &Modular<uint8_t>::write (std::ostream &os) const
	{
		return os << "Modular<uint8_t> mod " << uint32_t(_p);
	}

	std::istream &Modular<uint8_t>::read (std::istream &is)
	{
		is >> _p;
		return is;
	}

	std::ostream &Modular<uint8_t>::write (std::ostream &os, const Element &x) const
	{
		return os << int32_t(x);
	}

	std::istream &Modular<uint8_t>::read (std::istream &is, Element &x) const
	{
		int64_t tmp;
		is >> tmp;
		init(x,tmp);
		return is;
	}

} // namespace Givaro

#endif // __GIVARO_modular_uint8_INL

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
