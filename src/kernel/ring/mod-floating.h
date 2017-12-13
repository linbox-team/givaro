/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
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

/*! @file field/mod-floating.h
 * @ingroup field
 * @brief  representation of <code>Z/mZ</code> over \c float.
 */

#ifndef __GIVARO_mod_floating_H
#define __GIVARO_mod_floating_H

#include <float.h>

#include "givaro/givinteger.h"
#include "givaro/givcaster.h"
#include "givaro/givranditer.h"
#include "givaro/ring-interface.h"
#include "givaro/mod-general.h"
#include "givaro/mod-implem.h"

namespace Givaro
{

    template<typename _Storage_t, typename _Compute_t>
    class Mod<_Storage_t, _Compute_t, typename std::enable_if<std::is_floating_point<_Storage_t>::value>::type>: 
		public Mod_implem<_Storage_t, _Compute_t, typename make_unsigned_int<_Storage_t>::type>
    {
    public:

		using Storage_t = _Storage_t;
		using Compute_t = _Compute_t;
		using Residu_t = typename make_unsigned_int<_Storage_t>::type;

		using Element = Storage_t;
		using Self_t = Mod<_Storage_t, _Compute_t>;
		using Parent_t = Mod_implem<Storage_t, Compute_t, Residu_t>;

		// ----- Constructors
		using Parent_t::Mod_implem;

		using Parent_t::_p;
		using Parent_t::_pc;
		using Parent_t::_bitsizep;

        // ----- Exported Types and constantes
        //typedef Modular<float> Self_t;
        //typedef uint32_t Residu_t;
		//using Compute_t = float;
        //enum { size_rep = sizeof(Residu_t) };

        //// ----- Constantes
        //const Element zero = 0.f;
        //const Element one  = 1.f;
        //const Element mOne;

        //// ----- Constructors
        //Modular()
        //    : mOne(-1.f), _p(0.f), _lp(0)
        //{}

        //template<class T> Modular(const T& p)
        //    : mOne(Element(p) - 1.f), _p(Element(p)), _lp((Residu_t)p)
        //{
        //    assert(_p >= minCardinality());
        //    assert(_p <= maxCardinality());
        //}

        //Modular(const Self_t& F)
        //    : mOne(F.mOne), _p(F._p), _lp(F._lp)
        //{}

        //// ----- Accessors
        //inline Element minElement() const override { return zero; }
        //inline Element maxElement() const override { return mOne; }

        //// ----- Access to the modulus
        //inline Residu_t residu() const { return _lp; }
        //inline Residu_t size() const { return _lp; }
        //inline Residu_t characteristic() const { return _lp; }
        //inline float fcharacteristic() const { return _p; }
        //template<class T> inline T& characteristic(T& p) const { return p = _lp; }
        //inline Residu_t cardinality() const { return _lp; }
        //template<class T> inline T& cardinality(T& p) const { return p = _lp; }
		//static inline Residu_t maxCardinality() {
		//	// 2896 = max { p / 2*p^2 < 2^24 }
		//	return 2896;
		//}
        //static inline Residu_t minCardinality() { return 2; }

        //// ----- Checkers
        //inline bool isZero(const Element& a) const override { return a == zero; }
        //inline bool isOne (const Element& a) const override { return a == one; }
        //inline bool isMOne(const Element& a) const override { return a == mOne; }
        //inline bool isUnit(const Element& a) const override;
        //inline bool areEqual(const Element& a, const Element& b) const override { return a == b; }
        //inline size_t length(const Element a) const { return size_rep; }

        //// ----- Ring-wise operators
        //inline bool operator==(const Self_t& F) const { return _p == F._p; }
        //inline bool operator!=(const Self_t& F) const { return _p != F._p; }
        //inline Self_t& operator=(const Self_t& F)
        //{
        //    F.assign(const_cast<Element&>(one),  F.one);
        //    F.assign(const_cast<Element&>(zero), F.zero);
        //    F.assign(const_cast<Element&>(mOne), F.mOne);
        //    _p = F._p;
        //    _lp= F._lp;
        //    return *this;
        //}

        // ----- Initialisation
        Element& init (Element& x) const;
        
		__GIVARO_CONDITIONAL_TEMPLATE(Source, 
			std::is_same<Source,double>::value && 
			std::is_same<Element, float>::value)
		Element& init (Element& x, const Source y) const;
        
		__GIVARO_CONDITIONAL_TEMPLATE(Source, 
			std::is_integral<Source>::value && std::is_signed<Source>::value && 
			sizeof(Source) >= sizeof(Storage_t))
		Element& init (Element& x, const Source y) const;
        
		__GIVARO_CONDITIONAL_TEMPLATE(Source, 
			std::is_integral<Source>::value && std::is_unsigned<Source>::value && 
			sizeof(Source) >= sizeof(Storage_t))
		Element& init (Element& x, const Source y) const;
        
		__GIVARO_CONDITIONAL_TEMPLATE(Source, std::is_same<Source, Integer&>::value)
		Element& init (Element& x, const Source y) const;
		
		__GIVARO_CONDITIONAL_TEMPLATE(Source, 
			!(std::is_integral<Source>::value && sizeof(Source) >= sizeof(Storage_t)) &&
			!(std::is_same<Source, double>::value && std::is_same<Element, float>::value) &&
			!std::is_same<Source, Integer&>::value)
		Element& init(Element& r, const Source a) const
        { r = Caster<Element>(a); return reduce(r); }

        //Element& assign (Element& x, const Element& y) const;

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
        template< class Random > Element& random(Random& g, Element& r) const
        { return init(r, g()); }
        template< class Random > Element& nonzerorandom(Random& g, Element& a) const
        { while (isZero(init(a, g())))
                ;
            return a; }

        // --- IO methods
        //std::istream& read (std::istream& s);
        //std::ostream& write(std::ostream& s) const;
        //std::istream& read (std::istream& s, Element& a) const;
        //std::ostream& write(std::ostream& s, const Element& a) const;

    //protected:
    //    float _p;
    //    Residu_t _lp;
    };

} // Givaro

#include "mod-floating.inl"

#endif // __GIVARO_mod_floating_H

