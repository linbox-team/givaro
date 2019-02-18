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

namespace Givaro {

#define MOD Modular<Storage_t, Compute_t, typename std::enable_if<std::is_floating_point<Storage_t>::value>::type>
#define TMPL template<typename Storage_t, typename Compute_t>
#define COND_TMPL(T, ...) \
    template<typename T, \
    typename std::enable_if<(__VA_ARGS__), int>::type*>

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
    inline typename MOD::Element&
    MOD::init(Element& r, const Integer& a) const
    {
        r = Caster<Element>(a % _p);
        if (r < 0) r += _pc;
        return r;
    }

    // ------------
    // ----- Reduce

    TMPL
    inline typename MOD::Element&  MOD::reduce (Element& x) const
    {
        x = std::fmod(x, Caster<Element>(_pc));
        if (x < 0) x += Caster<Element>(_pc);
        return x;
    }

    TMPL
    inline typename MOD::Element&  MOD::reduce (Element& x, const Element& y) const
    {
        x = std::fmod(y, Caster<Element>(_pc));
        if (x < 0.f) x += Caster<Element>(_pc);
        return x;
    }

    // ------------------------
    // ----- Classic arithmetic

    TMPL
    inline typename MOD::Element & MOD::add
    (Element &x, const Element &y, const Element &z) const
    {
        Compute_t tmp = Caster<Compute_t>(y) + Caster<Compute_t>(z);
        return x = Caster<Element>(tmp<_pc?tmp:tmp-_pc);
    }

    TMPL
    inline typename MOD::Element & MOD::sub
    (Element &x, const Element &y, const Element &z) const
    {
        return x = (y>=z)? y-z:(Caster<Element>(_pc)-z)+y;
    }

    TMPL
    inline typename MOD::Element & MOD::mul
    (Element &x, const Element &y, const Element &z) const
    {
        return x = Caster<Element>(std::fmod(Caster<Compute_t>(y)*Caster<Compute_t>(z), _pc));
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
        return x = (y==Caster<Element>(0)?Caster<Element>(0):Caster<Element>(_pc)-y);
    }

    TMPL
    inline typename MOD::Element & MOD::inv
    (Element &x, const Element &y) const
    {
        invext(x,y,Caster<Element>(_pc));
        return (x<Caster<Element>(0) ? x+=Caster<Element>(_pc) : x);
    }

    TMPL
    inline typename MOD::Element & MOD::addin
    (Element &x, const Element &y) const
    {
        Compute_t tmp = Caster<Compute_t>(x) + Caster<Compute_t>(y);
        return x = Caster<Element>(tmp < _pc? tmp : tmp - _pc);
    }

    TMPL
    inline typename MOD::Element & MOD::subin
    (Element &x, const Element &y) const
    {
        if (x<y) x += (Caster<Element>(_pc)-y);
        else     x -= y;
        return x;
    }

    TMPL
    inline typename MOD::Element & MOD::mulin
    (Element &x, const Element &y) const
    {
        return x = Caster<Element>(std::fmod(Caster<Compute_t>(x)*Caster<Compute_t>(y), _pc));
    }

    TMPL
    inline typename MOD::Element & MOD::divin
    (Element &x, const Element &y) const
    {
        Element iy;
        return mulin(x, inv(iy, y));
    }

    TMPL
    inline typename MOD::Element & MOD::negin
    (Element &x) const
    {
        return x = (x==Caster<Element>(0)?Caster<Element>(0):Caster<Element>(_pc)-x);
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
        return r = Caster<Element>(std::fmod(Caster<Compute_t>(a)*Caster<Compute_t>(x)+Caster<Compute_t>(y), _pc));
    }

    TMPL
    inline typename MOD::Element & MOD::axpyin
    (Element &r, const Element &a, const Element &x) const
    {
        return r = Caster<Element>(std::fmod(Caster<Compute_t>(a)*Caster<Compute_t>(x)+Caster<Compute_t>(r), _pc));
    }

    // -- axmy: r <- a * x - y
    TMPL
    inline typename MOD::Element & MOD::axmy
    (Element& r, const Element &a, const Element &x, const Element &y) const
    {
        return r = Caster<Element>(std::fmod(Caster<Compute_t>(a)*Caster<Compute_t>(x) + (_pc-Caster<Compute_t>(y)), _pc));
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
        r = Caster<Element>(std::fmod(Caster<Compute_t>(a) * Caster<Compute_t>(x)
                                      + (_pc - Caster<Compute_t>(y)), _pc));
        return negin(r);
    }

    TMPL
    inline typename MOD::Element&  MOD::maxpyin
    (Element& r, const Element& a, const Element& x) const
    {
        Compute_t tmp = Caster<Compute_t>(a) * Caster<Compute_t>(x) + (_pc - Caster<Compute_t>(r));
        r = (tmp < _pc) ? Caster<Element>(tmp) : Caster<Element>(std::fmod(tmp, _pc));
        return negin(r);
    }

#undef MOD
#undef TMPL
#undef COND_TMPL

} // namespace Givaro

#endif // __GIVARO_modular_float_INL

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
