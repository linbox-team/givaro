// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: Brice Boyer (briceboyer) <boyer.brice@gmail.com>
//          A. Breust (taken from FFLAS-FFPACK)
// ==========================================================================

#ifndef __GIVARO_modular_double_INL
#define __GIVARO_modular_double_INL

#include "givaro/givcaster.h"
#include "givaro/modular-defines.h"

namespace Givaro {

	// --------------------
	// ----- Initialisation

	inline Modular<double>::Element& Modular<double>::init (Element& x) const
	{
		return x = zero;
	}

	inline Modular<double>::Element& Modular<double>::init (Element& x, const int64_t y) const
	{
		x = static_cast<Element>(std::abs(y) % _lp);
		if (y < 0) negin(x);
		return x;
	}

	inline Modular<double>::Element& Modular<double>::init (Element& x, const uint64_t y) const
	{
	    return x = static_cast<Element>(y % (uint64_t)(_lp));
	}

	inline Modular<double>::Element& Modular<double>::init (Element& x, const Integer& y) const
	{
	    x = static_cast<Element>(y % _lp);
            if (x < 0) x += _p;
            return x;
	}

	inline Modular<double>::Element& Modular<double>::assign (Element& x, const Element& y) const
	{
	    return x = y;
	}

	// ------------------------
	// ----- Convert and reduce

	inline Modular<double>::Element& Modular<double>::reduce (Element& x) const
	{
		x = fmod (x, _p);
		if (x < 0.0) x += _p;
		return x;
	}

	inline Modular<double>::Element& Modular<double>::reduce (Element& x, const Element& y) const
	{
		x = fmod (y, _p);
		if (x < 0.0) x += _p;
		return x;
	}

	// ------------------------
	// ----- Classic arithmetic

	inline Modular<double>::Element &Modular<double>::add
		(Element &x, const Element &y, const Element &z) const
	{
		__GIVARO_MODULAR_FLOATING_ADD(x,_p,y,z);
		return x;
	}

	inline Modular<double>::Element &Modular<double>::sub
		(Element &x, const Element &y, const Element &z) const
	{
		return __GIVARO_MODULAR_FLOATING_SUB(x,_p,y,z);
	}

	inline Modular<double>::Element &Modular<double>::mul
		(Element &x, const Element &y, const Element &z) const
	{
		return __GIVARO_MODULAR_FLOATING_MUL(x,_p,y,z);
	}

	inline Modular<double>::Element &Modular<double>::div
		(Element &x, const Element &y, const Element &z) const
	{
		return mulin(inv(x, z), y);
	}

	inline Modular<double>::Element &Modular<double>::neg
		(Element &x, const Element &y) const
	{
		return __GIVARO_MODULAR_FLOATING_NEG(x,_p,y);
	}

	inline Modular<double>::Element &Modular<double>::inv
		(Element &x, const Element &y) const
	{
		// The extended Euclidean algorithm
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
	}

	inline Modular<double>::Element &Modular<double>::addin
		(Element &x, const Element &y) const
	{
		__GIVARO_MODULAR_FLOATING_ADDIN(x,_p,y);
		return x;
	}

	inline Modular<double>::Element &Modular<double>::subin
		(Element &x, const Element &y) const
	{
		__GIVARO_MODULAR_FLOATING_SUBIN(x,_p,y);
		return x;
	}

	inline Modular<double>::Element &Modular<double>::mulin
		(Element &x, const Element &y) const
	{
		return __GIVARO_MODULAR_FLOATING_MULIN(x,_p,y);
	}

	inline Modular<double>::Element &Modular<double>::divin
		(Element &x, const Element &y) const
	{
		Modular<double>::Element iy;
		return mulin(x, inv(iy, y));
	}

	inline Modular<double>::Element &Modular<double>::negin
		(Element &x) const
	{
		return __GIVARO_MODULAR_FLOATING_NEGIN(x,_p);
	}

	inline Modular<double>::Element &Modular<double>::invin
		(Element &x) const
	{
		return inv(x, x);
	}

	// -- axpy: r <- a * x + y
	inline Modular<double>::Element &Modular<double>::axpy
		(Element &r, const Element &a, const Element &x, const Element &y) const
	{
		__GIVARO_MODULAR_FLOATING_MULADD(r, _p, a, x, y);
		return r;
	}

	inline Modular<double>::Element &Modular<double>::axpyin
		(Element &r, const Element &a, const Element &x) const
	{
		__GIVARO_MODULAR_FLOATING_MULADDIN(r, _p, a, x);
		return r;
	}

	// -- axmy: r <- a * x - y
	inline Modular<double>::Element &Modular<double>::axmy
		(Element& r, const Element &a, const Element &x, const Element &y) const
	{
		__GIVARO_MODULAR_FLOATING_MULSUB(r, _p, a, x, y);
		return r;
	}

	inline Modular<double>::Element &Modular<double>::axmyin
		(Element& r, const Element &a, const Element &x) const
	{
		maxpyin(r,a,x);
		return negin(r);
	}

	// -- maxpy:   r <- y - a * x
	inline Modular<double>::Element& Modular<double>::maxpy
		(Element& r, const Element& a, const Element& x, const Element& y) const
	{
		r = y;
		__GIVARO_MODULAR_FLOATING_SUBMULIN(r, _p, a, x);
		return r;
	}

	inline Modular<double>::Element& Modular<double>::maxpyin
		(Element& r, const Element& a, const Element& x) const
	{
		__GIVARO_MODULAR_FLOATING_SUBMULIN(r, _p, a, x);
		return r;
	}

	// ----------------
	// ----- IO methods

    inline
	std::ostream &Modular<double>::write (std::ostream &os) const
	{
		return os << "Modular<double> mod " << _lp;
	}

    inline
	std::istream &Modular<double>::read (std::istream &is)
	{
		is >> _p;
		return is;
	}

    inline
	std::ostream &Modular<double>::write (std::ostream &os, const Element &x) const
	{
            return os << x;
	}

    inline
	std::istream &Modular<double>::read (std::istream &is, Element &x) const
	{
		int64_t tmp;
		is >> tmp;
		init(x,tmp);
		return is;
	}

} // namespace Givaro

#endif // __GIVARO_modular_double_INL

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
