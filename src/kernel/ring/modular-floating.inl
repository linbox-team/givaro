// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: Brice Boyer (briceboyer) <boyer.brice@gmail.com>
//          A. Breust (taken from FFLAS-FFPACK)
// ==========================================================================

#ifndef __GIVARO_modular_float_INL
#define __GIVARO_modular_float_INL

#include "givaro/modular-defines.h"

namespace Givaro {

#define MOD Modular<Storage_t, Compute_t, typename std::enable_if<std::is_floating_point<Storage_t>::value>::type>
#define TMPL template<typename Storage_t, typename Compute_t> 
#define COND_TMPL(T, ...) \
	template<typename T, \
		typename std::enable_if<__VA_ARGS__, int>::type*>

    // --------------------
    // ----- Initialisation

    TMPL
    inline typename MOD::Element&
    MOD::init(Element& a) const
    {
        return a = this->zero;
    }

    TMPL
    COND_TMPL(Source, 
	std::is_same<Source, double>::value && std::is_same<Storage_t, float>::value)
    inline typename MOD::Element&
    MOD::init(Element& r, const Source a) const
    {
        r = Caster<Element>(std::fmod(a, _pc));
        if (r < 0.f) r += _pc;
        return r;
    }

    TMPL
    COND_TMPL(Source,
	std::is_integral<Source>::value && std::is_signed<Source>::value && sizeof(Source) >= sizeof(Storage_t))
    inline typename MOD::Element&
    MOD::init(Element& r, const Source a) const
    {
        r = Caster<Element>(std::abs(a) % Caster<Source>(_p));
        if (a < 0) negin(r);
        return r;
    }

    TMPL
    COND_TMPL(Source, 
	std::is_integral<Source>::value && std::is_unsigned<Source>::value && sizeof(Source) >= sizeof(Storage_t))
    inline typename MOD::Element&
    MOD::init(Element& r, const Source a) const
    {
        return r = Caster<Element>(a % Caster<Source>(_p));
    }

    TMPL
    COND_TMPL(Source, std::is_same<Source,  Integer&>::value)
    inline typename MOD::Element&
    MOD::init(Element& r, const Source a) const
    {
        r = Caster<Element>(a % _p);
        if (a < 0) negin(r);
        return r;
    }

    // ------------------------
    // ----- Convert and reduce

    TMPL
    inline typename MOD::Element&  MOD::reduce (Element& x) const
    {
        x = std::fmod(x, _pc);
        if (x < 0.f) x += _pc;
        return x;
    }

    TMPL
    inline typename MOD::Element&  MOD::reduce (Element& x, const Element& y) const
    {
        x = std::fmod(y, _pc);
        if (x < 0.f) x += _pc;
        return x;
    }

    // ------------------------
    // ----- Classic arithmetic

    TMPL
    inline typename MOD::Element & MOD::add
    (Element &x, const Element &y, const Element &z) const
    {
        __GIVARO_MODULAR_FLOATING_ADD(x,_pc,y,z);
        return x;
    }

    TMPL
    inline typename MOD::Element & MOD::sub
    (Element &x, const Element &y, const Element &z) const
    {
        return __GIVARO_MODULAR_FLOATING_SUB(x,_pc,y,z);
    }

    TMPL
    inline typename MOD::Element & MOD::mul
    (Element &x, const Element &y, const Element &z) const
    {
        return __GIVARO_MODULAR_FLOATING_MUL(x,_pc,y,z);
    }

    TMPL
    inline typename MOD::Element & MOD::div
    (Element &x, const Element &y, const Element &z) const
    {
        return mulin(inv(x, z), y);
    }

    TMPL
    inline typename MOD::Element & MOD::neg
    (Element &x, const Element &y) const
    {
        return __GIVARO_MODULAR_FLOATING_NEG(x,_pc,y);
    }

    TMPL
    inline typename MOD::Element & MOD::inv
    (Element &x, const Element &y) const
    {
        invext(x,y,Caster<Element>(_p));
        return (x<0 ? x+=Caster<Element>(_p) : x);
        return x;
	}

    //TMPL
    //inline bool MOD::isUnit(const Element& a) const 
    //{ 
    //    Element u,d; 
    //    invext(u,d,a,_pc); 
    //    return isOne(d) || isMOne(d); 
    //}

    TMPL
    inline typename MOD::Element & MOD::addin
    (Element &x, const Element &y) const
    {
        __GIVARO_MODULAR_FLOATING_ADDIN(x,_pc,y);
        return x;
    }

    TMPL
    inline typename MOD::Element & MOD::subin
    (Element &x, const Element &y) const
    {
        __GIVARO_MODULAR_FLOATING_SUBIN(x,_pc,y);
        return x;
    }

    TMPL
    inline typename MOD::Element & MOD::mulin
    (Element &x, const Element &y) const
    {
        return __GIVARO_MODULAR_FLOATING_MULIN(x,_pc,y);
    }

    TMPL
    inline typename MOD::Element & MOD::divin
    (Element &x, const Element &y) const
    {
        MOD::Element iy;
        return mulin(x, inv(iy, y));
    }

    TMPL
    inline typename MOD::Element & MOD::negin
    (Element &x) const
    {
        return __GIVARO_MODULAR_FLOATING_NEGIN(x,_pc);
    }

    TMPL
    inline typename MOD::Element & MOD::invin
    (Element &x) const
    {
        return inv(x, x);
    }

    // -- axpy: r <- a * x + y
    TMPL
    inline typename MOD::Element & MOD::axpy
    (Element &r, const Element &a, const Element &x, const Element &y) const
    {
        __GIVARO_MODULAR_FLOATING_MULADD(r, _pc, a, x, y);
        return r;
    }

    TMPL
    inline typename MOD::Element & MOD::axpyin
    (Element &r, const Element &a, const Element &x) const
    {
        __GIVARO_MODULAR_FLOATING_MULADDIN(r, _pc, a, x);
        return r;
    }

    // -- axmy: r <- a * x - y
    TMPL
    inline typename MOD::Element & MOD::axmy
    (Element& r, const Element &a, const Element &x, const Element &y) const
    {
        __GIVARO_MODULAR_FLOATING_MULSUB(r, _pc, a, x, y);
        return r;
    }

    TMPL
    inline typename MOD::Element & MOD::axmyin
    (Element& r, const Element &a, const Element &x) const
    {
        maxpyin(r,a,x);
        return negin(r);
    }

    // -- maxpy:   r <- y - a * x
    TMPL
    inline typename MOD::Element&  MOD::maxpy
    (Element& r, const Element& a, const Element& x, const Element& y) const
    {
        r = y;
        __GIVARO_MODULAR_FLOATING_SUBMULIN(r, _pc, a, x);
        return r;
    }

    TMPL
    inline typename MOD::Element&  MOD::maxpyin
    (Element& r, const Element& a, const Element& x) const
    {
        __GIVARO_MODULAR_FLOATING_SUBMULIN(r, _pc, a, x);
        return r;
    }

#undef MOD
#undef TMPL
#undef COND_TMPL

} // namespace Givaro

#endif // __GIVARO_modular_float_INL

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
