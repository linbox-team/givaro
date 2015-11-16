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

    // --------------------
    // ----- Initialisation

    inline Modular<float>::Element&
    Modular<float>::init(Element& a) const
    {
        return a = zero;
    }

    inline Modular<float>::Element&
    Modular<float>::init(Element& r, const double a) const
    {
        r = static_cast<float>(std::fmod(a, _p));
        if (r < 0.f) r += _p;
        return r;
    }

    inline Modular<float>::Element&
    Modular<float>::init(Element& r, const int32_t a) const
    {
        r = static_cast<Element>(std::abs(a) % _lp);
        if (a < 0) negin(r);
        return r;
    }

    inline Modular<float>::Element&
    Modular<float>::init(Element& r, const uint32_t a) const
    {
        return r = static_cast<Element>(a % uint32_t(_lp));
    }

    inline Modular<float>::Element&
    Modular<float>::init(Element& r, const int64_t a) const
    {
        r = static_cast<Element>(std::abs(a) % int64_t(_lp));
        if (a < 0) negin(r);
        return r;
    }

    inline Modular<float>::Element&
    Modular<float>::init(Element& r, const uint64_t a) const
    {
        return r = static_cast<Element>(a % uint64_t(_lp));
    }

    inline Modular<float>::Element&
    Modular<float>::init(Element& r, const Integer& a) const
    {
        r = static_cast<Element>(a % _lp);
        if (a < 0) negin(r);
        return r;
    }

    inline Modular<float>::Element&
    Modular<float>::assign (Element &x, const Element &y) const
    {
        return x = y;
    }

    // ------------------------
    // ----- Convert and reduce

    inline Modular<float>::Element& Modular<float>::reduce (Element& x) const
    {
        x = std::fmod(x, _p);
        if (x < 0.f) x += _p;
        return x;
    }

    inline Modular<float>::Element& Modular<float>::reduce (Element& x, const Element& y) const
    {
        x = std::fmod(y, _p);
        if (x < 0.f) x += _p;
        return x;
    }

    // ------------------------
    // ----- Classic arithmetic

    inline Modular<float>::Element &Modular<float>::add
    (Element &x, const Element &y, const Element &z) const
    {
        __GIVARO_MODULAR_FLOATING_ADD(x,_p,y,z);
        return x;
    }

    inline Modular<float>::Element &Modular<float>::sub
    (Element &x, const Element &y, const Element &z) const
    {
        return __GIVARO_MODULAR_FLOATING_SUB(x,_p,y,z);
    }

    inline Modular<float>::Element &Modular<float>::mul
    (Element &x, const Element &y, const Element &z) const
    {
        return __GIVARO_MODULAR_FLOATING_MUL(x,_p,y,z);
    }

    inline Modular<float>::Element &Modular<float>::div
    (Element &x, const Element &y, const Element &z) const
    {
        return mulin(inv(x, z), y);
    }

    inline Modular<float>::Element &Modular<float>::neg
    (Element &x, const Element &y) const
    {
        return __GIVARO_MODULAR_FLOATING_NEG(x,_p,y);
    }

    inline Modular<float>::Element &Modular<float>::inv
    (Element &x, const Element &y) const
    {
        // The extended Euclidean algorithm
        int32_t x_int, y_int, tx, ty;
        x_int = int32_t(_lp);
        y_int = int32_t(y);
        tx = 0;
        ty = 1;

        while (y_int != 0) {
            // always: gcd (modulus,residue) = gcd (x_int,y_int)
            //         sx*modulus + tx*residue = x_int
            //         sy*modulus + ty*residue = y_int
            int32_t q = x_int / y_int; // integer quotient
            int32_t temp = y_int;  y_int  = x_int  - q * y_int;
            x_int  = temp;
            temp = ty; ty = tx - q * ty;
            tx = temp;
        }

        if (tx < 0) tx += int32_t(_p);

        // now x_int = gcd (modulus,residue)
        return x = Element(tx);
    }

    inline Modular<float>::Element &Modular<float>::addin
    (Element &x, const Element &y) const
    {
        __GIVARO_MODULAR_FLOATING_ADDIN(x,_p,y);
        return x;
    }

    inline Modular<float>::Element &Modular<float>::subin
    (Element &x, const Element &y) const
    {
        __GIVARO_MODULAR_FLOATING_SUBIN(x,_p,y);
        return x;
    }

    inline Modular<float>::Element &Modular<float>::mulin
    (Element &x, const Element &y) const
    {
        return __GIVARO_MODULAR_FLOATING_MULIN(x,_p,y);
    }

    inline Modular<float>::Element &Modular<float>::divin
    (Element &x, const Element &y) const
    {
        Modular<float>::Element iy;
        return mulin(x, inv(iy, y));
    }

    inline Modular<float>::Element &Modular<float>::negin
    (Element &x) const
    {
        return __GIVARO_MODULAR_FLOATING_NEGIN(x,_p);
    }

    inline Modular<float>::Element &Modular<float>::invin
    (Element &x) const
    {
        return inv(x, x);
    }

    // -- axpy: r <- a * x + y
    inline Modular<float>::Element &Modular<float>::axpy
    (Element &r, const Element &a, const Element &x, const Element &y) const
    {
        __GIVARO_MODULAR_FLOATING_MULADD(r, _p, a, x, y);
        return r;
    }

    inline Modular<float>::Element &Modular<float>::axpyin
    (Element &r, const Element &a, const Element &x) const
    {
        __GIVARO_MODULAR_FLOATING_MULADDIN(r, _p, a, x);
        return r;
    }

    // -- axmy: r <- a * x - y
    inline Modular<float>::Element &Modular<float>::axmy
    (Element& r, const Element &a, const Element &x, const Element &y) const
    {
        __GIVARO_MODULAR_FLOATING_MULSUB(r, _p, a, x, y);
        return r;
    }

    inline Modular<float>::Element &Modular<float>::axmyin
    (Element& r, const Element &a, const Element &x) const
    {
        maxpyin(r,a,x);
        return negin(r);
    }

    // -- maxpy:   r <- y - a * x
    inline Modular<float>::Element& Modular<float>::maxpy
    (Element& r, const Element& a, const Element& x, const Element& y) const
    {
        r = y;
        __GIVARO_MODULAR_FLOATING_SUBMULIN(r, _p, a, x);
        return r;
    }

    inline Modular<float>::Element& Modular<float>::maxpyin
    (Element& r, const Element& a, const Element& x) const
    {
        __GIVARO_MODULAR_FLOATING_SUBMULIN(r, _p, a, x);
        return r;
    }

    // ----------------
    // ----- IO methods

    inline
    std::ostream &Modular<float>::write (std::ostream &os) const
    {
        return os << "Modular<float> mod " << _lp;
    }

    inline
    std::istream &Modular<float>::read (std::istream &is)
    {
        is >> _p;
        return is;
    }

    inline
    std::ostream &Modular<float>::write (std::ostream &os, const Element &x) const
    {
        return os << static_cast<uint64_t>(x);
    }

    inline
    std::istream &Modular<float>::read (std::istream &is, Element &x) const
    {
        int64_t tmp;
        is >> tmp;
        init(x,tmp);
        return is;
    }

} // namespace Givaro

#endif // __GIVARO_modular_float_INL

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
