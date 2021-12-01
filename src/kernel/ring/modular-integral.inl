// ==========================================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: Bruno Grenet, from existing files
// ==========================================================================
//

#ifndef __GIVARO_modular_integral_INL
#define __GIVARO_modular_integral_INL

namespace Givaro {

#define MOD Modular<Storage_t, Compute_t, \
    typename std::enable_if<std::is_integral<Storage_t>::value && std::is_integral<Compute_t>::value \
    && (sizeof(Storage_t) == sizeof(Compute_t) || 2*sizeof(Storage_t) == sizeof(Compute_t))>::type>
#define TMPL template<typename Storage_t, typename Compute_t>
#define COND_TMPL(T, ...) \
    template<typename T, \
    typename std::enable_if<(__VA_ARGS__), int>::type*>

    // ----- Initialisation

    TMPL
    inline typename MOD::Element&
    MOD::init (Element& x) const
    {
        return x = this->zero;
    }

    TMPL
    COND_TMPL(Source, IS_UINT(Source) && (sizeof(Source) > sizeof(Storage_t)))
    inline typename MOD::Element&
    MOD::init (Element& x, const Source y) const
    {
        x = Caster<Element>(y % Source(_p));
        return x;
    }

    TMPL
    COND_TMPL(Source, IS_SINT(Source) && (sizeof(Source) > sizeof(Storage_t)))
    inline typename MOD::Element&
    MOD::init (Element& x, const Source y) const
    {
        x = Caster<Element>((y<0 ? -y : y) % Source(_p));
        return (y < 0 ? negin(x) : x);
    }

    TMPL
    COND_TMPL(Source, IS_FLOAT(Source) && (sizeof(Source) >= sizeof(Storage_t)) && IS_SINT(Storage_t))
    inline typename MOD::Element&
    MOD::init (Element& x, const Source y) const
    {
        x = Caster<Element>(std::fmod(y, Source(_p)));
        if (x < Source(0.0)) x = Caster<Element>(x + _p);
        return x;
    }

    TMPL
    COND_TMPL(Source, IS_FLOAT(Source) && sizeof(Source) >= sizeof(Storage_t) && IS_UINT(Storage_t))
    inline typename MOD::Element&
    MOD::init (Element& x, const Source y) const
    {
        x = Caster<Element>(std::fmod((y < 0.0 ? -y : y), Source(_p)));
        return ( (y < 0.0) ? negin(x) : x);
    }

    TMPL
    inline typename MOD::Element&
    MOD::init (Element& x, const Integer& y) const
    {
        x = Caster<Element>(y % _p);
        if (x < 0) x = Caster<Element>(x + _p);
        return x;
    }

    TMPL
    COND_TMPL(Source, IS_UINT(Storage_t)
              &&!(IS_INT(Source) && (sizeof(Source) > sizeof(Storage_t)))
              &&!(IS_FLOAT(Source) && (sizeof(Source) >= sizeof(Storage_t))))
    inline typename MOD::Element&
    MOD::init (Element& x, const Source& y) const
    {
        reduce(x, Caster<Element>((y < 0)? -y : y));
        if (y < 0) negin(x);
        return x;
    }

    TMPL
    COND_TMPL(Source, IS_SINT(Storage_t)
              &&!(IS_INT(Source) && (sizeof(Source) > sizeof(Storage_t)))
              &&!(IS_FLOAT(Source) && (sizeof(Source) >= sizeof(Storage_t))))
    inline typename MOD::Element&
    MOD::init (Element& x, const Source& y) const
    {
        return reduce(Caster<Element>(x,y));
    }


    // ----- Reduce

    template<typename E, typename R, typename std::enable_if<IS_SINT(E), int>::type = 0>
    inline E& _reduce (E& x, const E& y, const R& p)
    {
        x = y % Caster<E>(p);
        return (x < 0 ? x = Caster<E>(x + p) : x);
    }

    template<typename E, typename R, typename std::enable_if<IS_UINT(E), int>::type = 0>
    inline E& _reduce (E& x, const E& y, const R& p)
    {
        return x = y % p;
    }

    template<typename E, typename R, typename std::enable_if<IS_SINT(E), int>::type = 0>
    inline E& _reduce (E& x, const R& p)
    {
        x %= Caster<E>(p);
        return (x < 0 ? x = Caster<E>(x + p) : x);
    }

    template<typename E, typename R, typename std::enable_if<IS_UINT(E), int>::type = 0>
    inline E& _reduce (E& x, const R& p)
    {
        return x %= p;
    }

    TMPL
    inline typename MOD::Element& MOD::reduce(Element& x, const Element& y) const
    {
        return x = _reduce<Element, Residu_t>(x, y, _p);
    }

    TMPL
    inline typename MOD::Element& MOD::reduce(Element& x) const
    {
        return x = _reduce<Element, Residu_t>(x, _p);
    }

    // ------------------------
    // ----- Classic arithmetic

    TMPL
    inline typename MOD::Element& MOD::mul
    (Element& r, const Element& a, const Element& b) const
    {
        return r = Caster<Element>(Caster<Compute_t>(a)*Caster<Compute_t>(b) % _pc);
    }

    TMPL
    inline typename MOD::Element& MOD::sub
    (Element& r, const Element& a, const Element& b) const
    {
        return r = (a < b) ? (Caster<Element>(_p) - b) + a : a - b;
    }

        // Add using overflows
    template<typename TElem, typename RElem,
             typename std::enable_if<! (IS_SINT(TElem)), int>::type = 0>
    inline TElem& GenericAdd(TElem& r, const TElem& a, const TElem& b, const RElem& _p)
    {
        r = a + b;
        return (r >= Caster<TElem>(_p) || r < a) ? r -= Caster<TElem>(_p) : r;
    }

        // Overflowing signed integrals is undefined
        // Thus clang++ 6.0.0, even though overflowing as expected,
        //   sometimes decides to not perform further operations...
        // The fix is to switch to unsigned values
        //   there overflowing is well defined.
        // Note, could swith to Compute_t if it is unsigned,
        //   but the unsigned type is probably sufficient.
    template<typename TElem, typename RElem,
             typename std::enable_if<IS_SINT(TElem), int>::type = 0>
    inline TElem& GenericAdd(TElem& r, const TElem& a, const TElem& b, const RElem& _p)
    {
        typename std::make_unsigned<TElem>::type rr(
            Caster<typename std::make_unsigned<TElem>::type>(a)
            +
            Caster<typename std::make_unsigned<TElem>::type>(b)
                                                    );
        return r = Caster<TElem>(
            rr >= Caster<typename std::make_unsigned<TElem>::type>(_p)
            || rr < Caster<typename std::make_unsigned<TElem>::type>(a) ?
            rr -= Caster<typename std::make_unsigned<TElem>::type>(_p)
            : rr);
    }

    TMPL
    inline typename MOD::Element& MOD::add
    (Element& r, const Element& a, const Element& b) const
    {
        return GenericAdd(r,a,b,_p);
    }

    TMPL
    inline typename MOD::Element& MOD::neg
    (Element& r, const Element& a) const
    {
        return r = (a == 0) ? Caster<Element>(0) : Caster<Element>(_p) - a;
    }

    TMPL
    inline typename MOD::Element& MOD::inv
    (Element& r, const Element& a) const
    {
        invext(r, a, Caster<Element>(_p));
        return (r < 0)? r += Caster<Element>(_p) : r;
    }

    TMPL
    inline typename MOD::Element& MOD::div
    (Element& r, const Element& a, const Element& b) const
    {
        return mulin(inv(r, b), a);
    }

    TMPL
    inline typename MOD::Element& MOD::mulin
    (Element& r, const Element& a) const
    {
        return r = Caster<Element>(Caster<Compute_t>(r)*Caster<Compute_t>(a) % _pc);
    }

    TMPL
    inline typename MOD::Element& MOD::divin
    (Element& r, const Element& a) const
    {
        Element ia;
        return mulin(r, inv(ia, a));
    }

        // Addin using overflows
    template<typename TElem, typename RElem,
             typename std::enable_if<! (IS_SINT(TElem)), int>::type = 0>
    inline TElem& GenericAddIN(TElem& r, const TElem& a, const RElem& _p)
    {
        r += a;
        return r = (r >= Caster<TElem>(_p) || r < a) ? r - Caster<TElem>(_p) : r;
    }

        // Addin using unsigned overflows, see comments for add
    template<typename TElem, typename RElem,
             typename std::enable_if<IS_SINT(TElem), int>::type = 0>
    inline TElem& GenericAddIN(TElem& r, const TElem& a, const RElem& _p)
    {
		typename std::make_unsigned<TElem>::type rr(r);
		rr += Caster<typename std::make_unsigned<TElem>::type>(a);
		return r = Caster<TElem>(
            rr >= Caster<typename std::make_unsigned<TElem>::type>(_p) ||
            rr < Caster<typename std::make_unsigned<TElem>::type>(a) ?
            rr -= Caster<typename std::make_unsigned<TElem>::type>(_p)
            : rr);
    }

    TMPL
    inline typename MOD::Element& MOD::addin
    (Element& r, const Element& a) const
    {
        return GenericAddIN(r,a,_p);
    }

    TMPL
    inline typename MOD::Element& MOD::subin
    (Element& r, const Element& a) const
    {
        return r = (r < a) ? (Caster<Element>(_p) - a) + r : r - a;
    }

    TMPL
    inline typename MOD::Element& MOD::negin
    (Element& r) const
    {
        return r = (r == 0) ? (Element)0 : Caster<Element>(_p) - r;
    }

    TMPL
    inline typename MOD::Element& MOD::invin
    (Element& r) const
    {
        return inv(r, r);
    }

    // -- axpy:   r <- a * x + y
    // -- axpyin: r <- a * x + r
    // -- axmy:   r <- a * x - y
    // -- axmyin: r <- a * x - r
    // -- maxpy:   r <- y - a * x
    // -- maxpyin: r <- r - a * x

    TMPL
    inline typename MOD::Element& MOD::axpy
    (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        return r = Caster<Element>((Caster<Compute_t>(a)*Caster<Compute_t>(b)
                                    + Caster<Compute_t>(c)) % _pc);
    }

    TMPL
    inline typename MOD::Element& MOD::axpyin
    (Element& r, const Element& a, const Element& b) const
    {
        return r = Caster<Element>((Caster<Compute_t>(a)*Caster<Compute_t>(b)
                                    + Caster<Compute_t>(r)) % _pc);
    }

    TMPL
    inline typename MOD::Element& MOD::maxpy
    (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        r = Caster<Element>((Caster<Compute_t>(a) * Caster<Compute_t>(b)
                             + (_pc - Caster<Compute_t>(c))) % _pc);
        return r = negin(r);
    }

    TMPL
    inline typename MOD::Element& MOD::axmy
    (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        return r = Caster<Element>((Caster<Compute_t>(a)*Caster<Compute_t>(b)
                                    + _pc-Caster<Compute_t>(c)) % _pc);
    }

    TMPL
    inline typename MOD::Element& MOD::maxpyin
    (Element& r, const Element& a, const Element& b) const
    {
        r = Caster<Element>((Caster<Compute_t>(a)*Caster<Compute_t>(b)
                             + _pc -Caster<Compute_t>(r)) % _pc);
        return r = negin(r);
    }

    TMPL
    inline typename MOD::Element& MOD::axmyin
    (Element& r, const Element& a, const Element& b) const
    {
        return r = Caster<Element>((Caster<Compute_t>(a)*Caster<Compute_t>(b)
                                    + _pc-Caster<Compute_t>(r)) % _pc);
    }

#undef MOD
#undef TMPL
#undef COND_TMPL

} // Givaro

#endif // __GIVARO_modular_integral_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
