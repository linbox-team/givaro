// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: A. Breust (adapted)
// ==========================================================================

#ifndef __GIVARO_modular_ruint_INL
#define __GIVARO_modular_ruint_INL

#include "modular-defines.h"

namespace Givaro
{
	// -------------
	// ----- Modular

    template<size_t K>
	inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Residu_t
	Modular<RecInt::ruint<K>, RecInt::ruint<K>>::getMaxModulus() { return RecInt::ruint<K>::getMaxModulus(); }

	// ------------------------
	// ----- Classic arithmetic

	template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::mul
		(Element& r, const Element& a, const Element& b) const
	{
		return __GIVARO_MODULAR_INTEGER_MUL(r,_p,a,b);
	}

	template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::sub
		(Element& r, const Element& a, const Element& b) const
	{
		return __GIVARO_MODULAR_INTEGER_SUB(r,_p,a,b);
	}

	template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::add
		(Element& r, const Element& a, const Element& b) const
	{
		__GIVARO_MODULAR_INTEGER_ADD(r,_p,a,b);
		return r;
	}

	template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::neg
		(Element& r, const Element& a) const
	{
		return __GIVARO_MODULAR_INTEGER_NEG(r,_p,a);
	}

	template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::inv
		(Element& r, const Element& a) const
	{
		return r = inv_mod(r, a, _p);
	}

	template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::div
		(Element& r, const Element& a, const Element& b) const
	{
		return mulin( inv(r,b), a );
	}

	template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::mulin
		(Element& r, const Element& a) const
	{
		return __GIVARO_MODULAR_INTEGER_MULIN(r,_p,a);
	}

	template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::divin
		(Element& r, const Element& a) const
	{
		Element ia;
		return mulin(r, inv(ia, a));
	}

	template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::addin
		(Element& r, const Element& a) const
	{
		__GIVARO_MODULAR_INTEGER_ADDIN(r,_p,a);
		return r;
	}

	template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::subin
		(Element& r, const Element& a) const
	{
		__GIVARO_MODULAR_INTEGER_SUBIN(r,_p,a);
		return r;
	}

	template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::negin
		(Element& r) const
	{
		return __GIVARO_MODULAR_INTEGER_NEGIN(r,_p);
	}

	template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::invin
		(Element& r) const
	{
		return inv(r, r);
	}

	template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::axpy
		(Element& r, const Element& a, const Element& b, const Element& c) const
	{
		return __GIVARO_MODULAR_INTEGER_MULADD(r,_p,a,b,c);
	}

	template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element&  Modular<RecInt::ruint<K>, RecInt::ruint<K>>::axpyin
		(Element& r, const Element& a, const Element& b) const
	{
		return __GIVARO_MODULAR_INTEGER_MULADDIN(r,_p,a,b);
	}

	template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::maxpy
		(Element& r, const Element& a, const Element& b, const Element& c) const
	{
		int16_t tmp;
		__GIVARO_MODULAR_INTEGER_MUL(tmp,_p,a,b);
		__GIVARO_MODULAR_INTEGER_SUB(r,_p,c,tmp);
		return r;
	}

	template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element&  Modular<RecInt::ruint<K>, RecInt::ruint<K>>::axmy
		(Element& r, const Element& a, const Element& b, const Element& c) const
	{
		return __GIVARO_MODULAR_INTEGER_MULSUB(r,_p,a,b,c);
	}

	template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element&  Modular<RecInt::ruint<K>, RecInt::ruint<K>>::maxpyin
		(Element& r, const Element& a, const Element& b) const
	{
		__GIVARO_MODULAR_INTEGER_SUBMULIN(r,_p,a,b);
		return r;
	}

	template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element&  Modular<RecInt::ruint<K>, RecInt::ruint<K>>::axmyin
		(Element& r, const Element& a, const Element& b) const
	{
		return __GIVARO_MODULAR_INTEGER_MULSUB(r,_p,a,b,r);
	}

	// --------------------
	// ----- Initialisation

	template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element&  Modular<RecInt::ruint<K>, RecInt::ruint<K>>::init
        (Element& r, const unsigned long a) const
	{
		return r = (Element)( a >= (unsigned long)_p ? a % (unsigned long)_p : a);
	}

	template<size_t K>
    inline  typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element&  Modular<RecInt::ruint<K>, RecInt::ruint<K>>::init ( Element& r, const long a ) const
	{
		int sign; long ua;
		if (a <0) {
			sign =-1;
			ua = -a;
		}
		else {
			ua = a;
			sign =1;
		}
		r = Element( (ua >=_p) ? ua % (uint16_t)_p : ua );
		if (r && (sign ==-1))
			r = (Element)(_p - r);
		return r;
	}

	template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element&  Modular<RecInt::ruint<K>, RecInt::ruint<K>>::init ( Element& r, const Integer& Residu ) const
	{
		Element tr;
		if (Residu <0) {
			// -a = b [p]
			// a = p-b [p]
			if ( Residu <= (Integer)(-_p) ) tr = Element( (-Residu) % (uint16_t)_p) ;
			else tr = Element(-Residu);
			if (tr)
				return r = Element((uint16_t)_p - (uint16_t)tr);
			else
				return r = zero;
		}
		else {
			if (Residu >= (Integer)_p ) tr =   Element(Residu % _p) ;
			else tr = Element(Residu);
			return r = (Element)tr;
		}
	}

	template<size_t K>
    inline  typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::init( Element& a, const int i) const
	{
		return init(a,(long)i);
	}
	template<size_t K>
    inline  typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::init( Element& a, const double i) const
	{
		return init(a,(long)i);
	}
	template<size_t K>
    inline  typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::init( Element& a, const float i) const
	{
		return init(a,(double)i);
	}
	template<size_t K>
    inline  typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::init( Element& a, const unsigned int i) const
	{
		return init(a,(unsigned long)i);
	}

	template<size_t K>
    inline  typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element&  Modular<RecInt::ruint<K>, RecInt::ruint<K>>::assign ( Element& r, const Element a ) const
	{  return r=a;
	}

	template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::init ( Element& r ) const
	{
		return r = zero;
	}

	// -- Input: (z, <_p>)
	template<size_t K>
	inline std::ostream& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::write (std::ostream& s) const
	{
		return s << "Modular<RecInt::ruint<K>, RecInt::ruint<K>> modulo " << residu();
	}

	template<size_t K>
    inline std::istream& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::read (std::istream& s, Element& a) const
	{
        Integer tmp;
		s >> tmp;
		init(a, tmp);
		return s;
	}

	template<size_t K>
    inline std::ostream& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::write (std::ostream& s, const Element a) const
	{
		return s << a;
	}
}

#endif
