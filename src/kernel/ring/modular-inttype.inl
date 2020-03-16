// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givzpzGen.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: JG Dumas
// $Id: givzpzGen.inl,v 1.11 2011-02-02 17:16:43 briceboyer Exp $
// ==========================================================================
#ifndef __GIVARO_zpz_gen_INL
#define __GIVARO_zpz_gen_INL
// Description:

namespace Givaro
{
    // ------------------------- Arithmetic functions

    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::mul (Element& r, const Element& a, const Element& b) const
    {
        r = a;
        r *= b;
        r %= _p;
        return r;
    }

    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::sub (Element& r, const Element& a, const Element& b) const
    {
        return r = (a>=b)? a-b : (_p-b) + a;
    }

    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::add (Element& r, const Element& a, const Element& b) const
    {
        r = a+b;
        if (r >= _p) r -= _p;
        return r;
    }

    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::neg (Element& r, const Element& a) const
    {
        return r = (a == zero ? zero : _p - a);
    }

    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::negin (Element& r) const
    {
        return r = (Modular<IntType, COMP, Enable>::isZero(r) ? zero : _p - r);
    }

    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::inv (Element& u1, const Element& a) const
    {
        u1=one;
        IntType r0(_p), r1(a);
        IntType q(r0/r1);

        r0 -= q * r1;
        if (r0 == zero) return u1;
        IntType u0 = q;

        q = r1/r0;
        r1 -= q * r0;

        while (r1 != zero) {
            u1 += q * u0;

            q = r0/r1;
            r0 -= q * r1;
            if (r0 == zero) return u1;
            u0 += q * u1;

            q = r1/r0;
            r1 -= q * r0;

        };

        return u1=_p-u0;
    }

    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::div (Element& r, const Element& a, const Element& b) const
    {
        typename Modular<IntType, COMP, Enable>::Element ib;
        inv(ib, b);
        r = a; r*=ib; r %= _p;
        return r;
    }

    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::mulin (Element& r, const Element& a) const
    {
        r *= a;
        r %= _p;
        return r;
    }

    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::divin (Element& r, const Element& a) const
    {
        typename Modular<IntType, COMP, Enable>::Element ia;
        inv(ia, a);
        return mulin(r, ia);
    }

    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::addin (Element& r, const Element& a) const
    {
        r += a;
        if (r >= _p) r-=_p;
        return r;
    }

    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::subin (Element& r, const Element& a) const
    {
        if (r<a) r += (_p - a);
        else     r -= a;
        return r;
    }


    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::invin (Element& r) const
    {
        typename Modular<IntType, COMP, Enable>::Element t = r;
        return Modular<IntType, COMP, Enable>::inv(r,t);
    }

    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::axpy (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        r = a;
        r *= b;
        r += c;
        r %= _p;
        return r;
    }

    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::axpyin (Element& r, const Element& a, const Element& b) const
    {
        typename Modular<IntType, COMP, Enable>::Element tmp = r;
        tmp += a*b;
        tmp %= _p;
        return r = (Element)tmp;
    }

    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::axmy (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        r = a*b;
        r += (_p - c);
        return r = (r < _p ? r : r%_p);
    }

    // r = c - a*b
    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::maxpy (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        r = c;
        return maxpyin(r, a, b);
    }

    // r -= a*b
    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::maxpyin (Element& r, const Element& a, const Element& b) const
    {
        axmyin(r,a,b);
        return negin(r);
    }
    // r = a*b - r
    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::axmyin (Element& r, const Element& a, const Element& b) const
    {
        r = _p - r;
        r += a*b;
        return r = (r<_p ? r : r % _p);
    }


} // namespace Givaro

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
