/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
// ==========================================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: Bruno Grenet, from existing files
// ==========================================================================
//

/*! @file ring/mod-integral.h
 * @ingroup ring
 * @brief  representation of <code>Z/mZ</code> over int types.
 */

#ifndef __GIVARO_mod_integral_H
#define __GIVARO_mod_integral_H

#include "givaro/givinteger.h"
#include "givaro/givcaster.h"
#include "givaro/givranditer.h"
#include "givaro/ring-interface.h"
#include "givaro/mod-implem.h"
//#include "modular-defines.h"

#define VALID_TYPES(S, C) \
	typename std::enable_if<std::is_integral<S>::value && std::is_integral<C>::value \
		&& (sizeof(S) == sizeof(C) || 2*sizeof(S) == sizeof(C))>::type 
#define IS_INT(T) std::is_integral<T>::value
#define IS_SINT(T) std::is_integral<T>::value && std::is_signed<T>::value
#define IS_UINT(T) std::is_integral<T>::value && std::is_unsigned<T>::value
#define IS_FLOAT(T) std::is_floating_point<T>::value
#define IS(S,T) std::is_same<S, T>::value

namespace Givaro {

	// -- Enabled for (Storage_t, Compute_t) either ((u)intN_t, (u)intN_t) or ((u)intN_t, (u)int2N_t)
	// -- Note that Compute_t is always converted to its unsigned version


	template<typename _Storage_t, typename _Compute_t>
	class Mod<_Storage_t, _Compute_t, VALID_TYPES(_Storage_t, _Compute_t)>: 
		public Mod_implem<_Storage_t, typename std::make_unsigned<_Compute_t>::type, typename std::make_unsigned<_Storage_t>::type>
	{
	public:
		
		using Storage_t = _Storage_t;
		using Residu_t = typename std::make_unsigned<_Storage_t>::type;
		using Compute_t = typename std::make_unsigned<_Compute_t>::type;

		using Element = Storage_t;
		using Self_t = Mod<Storage_t, _Compute_t>;
		using Parent_t = Mod_implem<Storage_t, Compute_t, Residu_t>;

        // ----- Constructors
		using Parent_t::Mod_implem;

		using Parent_t::_p;
		using Parent_t::_pc;
		//using Parent_t::_bitsizep;

    	//inline bool isUnit(const Element& a) const 
    	//{ 
    	//    Element u,d; 
    	//    invext(u,d,a,Element(_p)); 
    	//    return isOne(d) || isMOne(d); 
    	//}
		
		// ----- Initialisation

		Element& init (Element& x) const
		{
			return x = this->zero;
		}

		__GIVARO_CONDITIONAL_TEMPLATE(Source, IS_UINT(Source) && (sizeof(Source) > sizeof(Element)))
		Element& init (Element& x, const Source y) const
		{
			x = Caster<Element>(y % Source(_p));
			return x;
		}

		__GIVARO_CONDITIONAL_TEMPLATE(Source, IS_SINT(Source) && (sizeof(Source) > sizeof(Element)))
		Element& init (Element& x, const Source y) const
		{
			x = Caster<Element>((y<0 ? -y : y) % Source(_p));
			return (y < 0 ? negin(x) : x);
		}

		__GIVARO_CONDITIONAL_TEMPLATE(Source, IS_FLOAT(Source) && (sizeof(Source) >= sizeof(Element)) && IS_SINT(Element))
		Element& init (Element& x, const Source y) const
		{
			x = Caster<Element>(std::fmod(y, Source(_p)));
			if (x < Source(0.0)) x = Caster<Element>(x + _p);
			return x;
		}

		__GIVARO_CONDITIONAL_TEMPLATE(Source, IS_FLOAT(Source) && sizeof(Source) >= sizeof(Element) && IS_UINT(Element))
		Element& init (Element& x, const Source y) const
		{
			x = Caster<Element>(std::fmod((y < 0.0 ? -y : y), Source(_p)));
			if (x < Source(0.0)) x = Caster<Element>(x + _p);
			return x;
		}

		__GIVARO_CONDITIONAL_TEMPLATE(Source, IS(Source, Integer&))
		Element& init (Element& x, const Source y) const
		{
			x = Caster<Element>(y % _p);
			if (x < 0) x = Caster<Element>(x + _p);
			return x;
		}

		__GIVARO_CONDITIONAL_TEMPLATE(Source, IS_UINT(Element) 
						&&!(IS_INT(Source) && (sizeof(Source) > sizeof(Element)))
						&&!(IS_FLOAT(Source) && (sizeof(Source) >= sizeof(Element)))
						&&!IS(Source, Integer&))
		Element& init (Element& x, const Source y) const
		{
			reduce(x, Caster<Element>((y < 0)? -y : y));
			if (y < 0) negin(x);
			return x;
		}

		__GIVARO_CONDITIONAL_TEMPLATE(Source, IS_SINT(Element) 
						&&!(IS_INT(Source) && (sizeof(Source) > sizeof(Element)))
						&&!(IS_FLOAT(Source) && (sizeof(Source) >= sizeof(Element)))
						&&!IS(Source, Integer&))
		Element& init (Element& x, const Source y) const
		{
			return reduce(Caster<Element>(x,y)); 
		}


		// ----- Reduce

		//template<typename E=Element>
		//typename std::enable_if<std::is_signed<E>::value, E&>::type
		__GIVARO_CONDITIONAL_TEMPLATE(E = Element, IS_SINT(E))
		E& reduce (E& x, const E& y) const
		{
			x = y % Caster<E>(_p);
			return (x < 0 ? x = Caster<E>(x + _p) : x);
		}

		//template<typename E=Element>
		//typename std::enable_if<std::is_unsigned<E>::value, E&>::type
		__GIVARO_CONDITIONAL_TEMPLATE(E = Element, IS_UINT(E))
		E& reduce (E& x, const E& y) const
		{
			return x = y % _p;
		}

		//template<typename E=Element>
		//typename std::enable_if<std::is_signed<E>::value, E&>::type
		__GIVARO_CONDITIONAL_TEMPLATE(E = Element, IS_SINT(E))
		E& reduce (E& x) const
		{
			x %= Caster<E>(_p);
			return (x < 0 ? x = Caster<E>(x + _p) : x);
		}

		//template<typename E=Element>
		//typename std::enable_if<std::is_unsigned<E>::value, E&>::type
		__GIVARO_CONDITIONAL_TEMPLATE(E = Element, IS_UINT(E))
		E& reduce (E& x) const
		{
			return x %= _p;
		}

		// ------------------------
		// ----- Classic arithmetic
	
		inline Element&
		mul
		(Element& r, const Element& a, const Element& b) const
		{
			return r = Caster<Element>(Caster<Compute_t>(a)*Caster<Compute_t>(b) % Caster<Compute_t>(_p));
		}
	
		inline Element&
		sub
		(Element& r, const Element& a, const Element& b) const
		{
			return r = Caster<Element>((Caster<Compute_t>(a) >= Caster<Compute_t>(b)) ?
				Caster<Compute_t>(a)-Caster<Compute_t>(b) :
				Caster<Compute_t>(_p)-Caster<Compute_t>(b)+Caster<Compute_t>(a));
			//TODO: Caster necessary? define variables...
		}
	
		inline Element&
		add
		(Element& r, const Element& a, const Element& b) const
		{
    		Compute_t tmp = Caster<Compute_t>(a) + Caster<Compute_t>(b); 
    		return r = Caster<Element>((tmp < Caster<Compute_t>(_p)) ? tmp : tmp - Caster<Compute_t>(_p));
		}
	
		inline Element&
		neg
		(Element& r, const Element& a) const
		{
			return r = (a == 0) ? Caster<Element>(0) : Caster<Element>(Caster<Compute_t>(_p)-Caster<Compute_t>(a));
		}
	
		inline Element&
		inv
		(Element& r, const Element& a) const
		{
			invext(r, a, Caster<Element>(_p));
			return (r < 0)? r += _p : r;
		}
	
		inline Element&
		div
		(Element& r, const Element& a, const Element& b) const
		{
			return mulin(inv(r, b), a);
		}
	
		inline Element&
		mulin
		(Element& r, const Element& a) const
		{
			return r = Caster<Element>(Caster<Compute_t>(r)*Caster<Compute_t>(a) % Caster<Compute_t>(_p));
		}
	
		inline Element&
		divin
		(Element& r, const Element& a) const
		{
			Element ia;
			return mulin(r, inv(ia, a));
		}
	
		inline Element&
		addin
		(Element& r, const Element& a) const
		{
			Compute_t tmp = Caster<Compute_t>(r) + Caster<Compute_t>(a);
			return r = Caster<Element>((tmp < _p) ? tmp : tmp - _p);
			//TODO: No need to cast to Compute_t: detect overflow
		}
	
		inline Element&
		subin
		(Element& r, const Element& a) const
		{
			Compute_t rc = Caster<Compute_t>(r), ac = Caster<Compute_t>(a);
			return r = Caster<Element>((rc < ac) ? Caster<Compute_t>(_p) + rc - ac : rc - ac);
		}
	
		inline Element&
		negin
		(Element& r) const
		{
			return r = (r == 0) ? (Element)0 : Caster<Element>(Caster<Compute_t>(_p)-Caster<Compute_t>(r));
		}
	
		inline Element&
		invin
		(Element& r) const
		{
			return inv(r, r);
		}

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
		// -- axmy:   r <- a * x - y
		// -- axmyin: r <- a * x - r
		// -- maxpy:   r <- y - a * x
		// -- maxpyin: r <- r - a * x
	
		inline Element&
		axpy
		(Element& r, const Element& a, const Element& b, const Element& c) const
		{
			return r = Caster<Element>((Caster<Compute_t>(a)*Caster<Compute_t>(b) 
						+ Caster<Compute_t>(c)) % Caster<Compute_t>(_p));
		}
	
		inline Element&
		axpyin
		(Element& r, const Element& a, const Element& b) const
		{
			return r = Caster<Element>((Caster<Compute_t>(a)*Caster<Compute_t>(b) 
						+ Caster<Compute_t>(r)) % Caster<Compute_t>(_p));
		}
	
		inline Element&
		maxpy
		(Element& r, const Element& a, const Element& b, const Element& c) const
		{
			Element tmp;
			return r = sub(r, c, mul(tmp, a, b));
			// TODO: Use only one modulo?
		}
	
		inline Element&
		axmy
		(Element& r, const Element& a, const Element& b, const Element& c) const
		{
    		return r = Caster<Element>((Caster<Compute_t>(a)*Caster<Compute_t>(b)
    					+ Caster<Compute_t>(_p)-Caster<Compute_t>(c)) % Caster<Compute_t>(_p));
		}
	
		inline Element&
		maxpyin
		(Element& r, const Element& a, const Element& b) const
		{
    		r = Caster<Element>((Caster<Compute_t>(a)*Caster<Compute_t>(b) 
				+ Caster<Compute_t>(_p) -Caster<Compute_t>(r)) % Caster<Compute_t>(_p));
			return r = negin(r);
		}
	
		inline Element&
		axmyin
		(Element& r, const Element& a, const Element& b) const
		{
    		return r = Caster<Element>((Caster<Compute_t>(a)*Caster<Compute_t>(b)
    					+ Caster<Compute_t>(_p)-Caster<Compute_t>(r)) % Caster<Compute_t>(_p));
		}

		// ----- Random generators
		typedef ModularRandIter<Self_t> RandIter;
		typedef GeneralRingNonZeroRandIter<Self_t> NonZeroRandIter;
		template< class Random > Element& random(Random& g, Element& r) const
		{ return init(r, g()); }
		template< class Random > Element& nonzerorandom(Random& g, Element& a) const
		{ while (isZero(init(a, g())))
				;
			return a; }

	};
}

#undef VALID_TYPES
#undef IS_INT
#undef IS_SINT
#undef IS_UINT
#undef IS_FLOAT
#undef IS

#endif // __GIVARO_mod_integral_H
