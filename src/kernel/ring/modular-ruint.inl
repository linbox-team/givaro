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

namespace Givaro
{
#define TMPL template<typename Storage_t, typename Compute_t>
#define MOD Modular<Storage_t, Compute_t, \
    typename std::enable_if<is_same_ruint<Storage_t, Compute_t>::value \
    || is_smaller_ruint<Storage_t, Compute_t>::value \
    || is_same_rint<Storage_t, Compute_t>::value \
    || is_smaller_rint<Storage_t, Compute_t>::value>::type>

#define DIFF_RECINT \
    typename std::enable_if<!std::is_same<E,C>::value, int>::type = 0
#define SAME_RECINT \
    typename std::enable_if<std::is_same<E,C>::value, int>::type = 0

    template<typename E, typename C, DIFF_RECINT>
    E& _mul
    (E& r, const E& a, const E& b, const E& p)
    {
        C tmp;
        RecInt::lmul(tmp, a, b);
        RecInt::mod_n(r, tmp, p);
        return r;
    }

    template<typename E, typename C, SAME_RECINT>
    E& _mul
    (E& r, const E& a, const E& b, const E& p)
    {
        E tmp;
        RecInt::mul(tmp, a, b);
        RecInt::mod_n(r, tmp, p);
        return r;
    }

    TMPL
    inline
    typename MOD::Element& MOD::mul
    (Element& r, const Element& a, const Element& b) const
    {
        return _mul<Element, Compute_t>(r, a, b, _p);
    }

    TMPL
    inline typename MOD::Element& MOD::sub
    (Element& r, const Element& a, const Element& b) const
    {
        if (a < b) {
            RecInt::sub(r, _p, b);
            RecInt::add(r, a);
        }
        else RecInt::sub(r, a, b);
        return r;
    }

    TMPL
    inline typename MOD::Element& MOD::add
    (Element& r, const Element& a, const Element& b) const
    {
        RecInt::add(r, a, b);
        if (r >= _p) RecInt::sub(r, _p);
        return r;
    }

    TMPL
    inline typename MOD::Element& MOD::neg
    (Element& r, const Element& a) const
    {
        if (a == 0) RecInt::reset(r);
        else        RecInt::sub(r, _p, a);
        return r;
    }

    TMPL
    inline typename MOD::Element& MOD::inv
    (Element& r, const Element& a) const
    {
        return inv_mod(r, a, _p);
    }

    TMPL
    inline typename MOD::Element& MOD::div
    (Element& r, const Element& a, const Element& b) const
    {
        return mulin( inv(r,b), a );
    }

    template<typename E, typename C>
    E& _mulin(E& r, const E& a, const E& p, const C& pc)
    {
        C tmp;
        RecInt::lmul(tmp, r, a);
        RecInt::mod_n(r, tmp, p);
        return r;
    }
    template<typename E>
    E& _mulin(E& r, const E& a, const E& p, const E& pc)
    {
        RecInt::mod_n(RecInt::mul(r, a), p);
        return r;
    }

    TMPL
    inline
    typename MOD::Element& MOD::mulin
    (MOD::Element& r, const MOD::Element& a) const
    {
        return _mulin(r, a, _p, _pc);
    }

    TMPL
    inline typename MOD::Element& MOD::divin
    (Element& r, const Element& a) const
    {
        Element ia;
        return mulin(r, inv(ia, a));
    }

    TMPL
    inline typename MOD::Element& MOD::addin
    (Element& r, const Element& a) const
    {
        RecInt::add(r, a);
        if (r >= _p) RecInt::sub(r, _p);
        return r;
    }

    TMPL
    inline typename MOD::Element& MOD::subin
    (Element& r, const Element& a) const
    {
        if (r < a) RecInt::add(r, _p - a);
        else       RecInt::sub(r, a);
        return r;
    }

    TMPL
    inline typename MOD::Element& MOD::negin
    (Element& r) const
    {
        if (r == 0) RecInt::reset(r);
        else        RecInt::sub(r, _p, r);
        return r;
    }

    TMPL
    inline typename MOD::Element& MOD::invin
    (Element& r) const
    {
        return inv(r, r);
    }

    template<typename E, typename C, DIFF_RECINT>
    inline E& _axpy (E& r, const E& a, const E& b, const E& c, const E& p)
    {
        C tmp;
        RecInt::lmul(tmp, a, b);
        RecInt::mod_n(r, tmp, p);
        RecInt::add(r, c);
        if (r >= p) RecInt::sub(r, p);
        return r;
    }
    template<typename E, typename C, SAME_RECINT>
    inline E& _axpy (E& r, const E& a, const E& b, const E& c, const E& p)
    {
        RecInt::copy(r, c);
        RecInt::addmul(r, a, b);
        RecInt::mod_n(r, p);
        return r;
    }

    TMPL
    inline typename MOD::Element& MOD::axpy
    (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        return _axpy<Element, Compute_t>(r, a, b, c, _p);
    }

    template<typename E, typename C, DIFF_RECINT>
    inline E& _axpyin (E& r, const E& a, const E& b, const E& p)
    {
        E tmp = r;
        return r = _axpy<E, C>(r, a, b, tmp, p);
    }
    template<typename E, typename C, SAME_RECINT>
    inline E& _axpyin (E& r, const E& a, const E& b, const E& p)
    {
        RecInt::addmul(r, a, b);
        RecInt::mod_n(r, p);
        return r;
    }
    TMPL
    inline typename MOD::Element&  MOD::axpyin
    (Element& r, const Element& a, const Element& b) const
    {
        return _axpyin<Element, Compute_t>(r, a, b, _p);
    }

    TMPL
    inline typename MOD::Element& MOD::maxpy
    (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        _mul<Element, Compute_t>(r, a, b, _p);
        sub(r, c, r);
        return r;
    }

    TMPL
    inline typename MOD::Element&  MOD::axmy
    (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        _mul<Element, Compute_t>(r, a, b, _p);
        subin(r, c);
        return r;
    }

    template<typename E, typename C, DIFF_RECINT>
    E& _maxpyin (E& r, const E& a, const E& b, const E& p)
    {
        E tmp;
        _mul<E, C>(tmp, a, b, p);
        if (r < tmp)
        {
            RecInt::sub(tmp, p, tmp);
            RecInt::add(r, tmp);
        }
        else RecInt::sub(r, tmp);
        return r;
    }
    template<typename E, typename C, SAME_RECINT>
    E& _maxpyin (E& r, const E& a, const E& b, const E& p)
    {
        if (r == 0) RecInt::reset(r);
        else RecInt::sub(r, p, r);
        RecInt::addmul(r, a, b);
        RecInt::mod_n(r, p);
        if (r == 0) RecInt::reset(r);
        else RecInt::sub(r, p, r);
        return r;
    }
    TMPL
    inline typename MOD::Element&  MOD::maxpyin
    (Element& r, const Element& a, const Element& b) const
    {
        return _maxpyin<Element, Compute_t> (r, a, b, _p);
    }

    TMPL
    inline typename MOD::Element&  MOD::axmyin
    (Element& r, const Element& a, const Element& b) const
    {
        Element rc(r);
        axmy(r, a, b, rc);
        return r;
    }

#undef TMPL
#undef MOD
#undef SAME_RECINT
#undef DIFF_RECINT

}


#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
