// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: Brice Boyer (briceboyer) <boyer.brice@gmail.com>
//          A. Breust (taken from FFLAS-FFPACK)
// ==========================================================================

#ifndef __GIVARO_modular_uint8_INL
#define __GIVARO_modular_uint8_INL

#include <cmath> // fmod

#include "givaro/modular-defines.h"

namespace Givaro {

    // -------------
    // ----- Modular

    template<>
    inline Modular<uint8_t, int8_t>::Residu_t
    Modular<uint8_t, int8_t>::maxCardinality() { return 15u; } // 2^4 - 1

    template<>
    inline Modular<uint8_t, uint8_t>::Residu_t
    Modular<uint8_t, uint8_t>::maxCardinality() { return 15u; }

    template<>
    inline Modular<uint8_t, uint16_t>::Residu_t
    Modular<uint8_t, uint16_t>::maxCardinality() { return 127u; } // 2^7 - 1

    template<>
    inline Modular<uint8_t, int16_t>::Residu_t
    Modular<uint8_t, int16_t>::maxCardinality() { return 127u; }

    // --------------------
    // ----- Initialisation

    template<typename COMP> inline typename Modular<uint8_t, COMP>::Element&
    Modular<uint8_t, COMP>::init (Element& r, const float a) const
    {
	r = static_cast<Element>(std::fmod((a < 0.f)? -a : a, float(_p)));
	return r = (a < 0.f)? negin(r) : r;
    }

    template<typename COMP> inline typename Modular<uint8_t, COMP>::Element&
    Modular<uint8_t, COMP>::init (Element& r, const double a) const
    {
	r = static_cast<Element>(std::fmod((a < 0.0)? -a : a, double(_p)));
	return r = (a < 0.0)? negin(r) : r;
    }

    template<typename COMP> inline typename Modular<uint8_t, COMP>::Element&
    Modular<uint8_t, COMP>::init (Element& r, const int16_t a) const
    {
	r = static_cast<Element>(((a < 0)? -a : a) % int16_t(_p));
        if (a < 0) negin(r);
	return r;
    }

    template<typename COMP> inline typename Modular<uint8_t, COMP>::Element&
    Modular<uint8_t, COMP>::init (Element& r, const uint16_t a) const
    {
	r = static_cast<Element>(a % uint16_t(_p));
	return r;
    }

    template<typename COMP> inline typename Modular<uint8_t, COMP>::Element&
    Modular<uint8_t, COMP>::init (Element& r, const int32_t a) const
    {
	r = static_cast<Element>(std::abs(a) % int32_t(_p));
        if (a < 0) negin(r);
	return r;
    }

    template<typename COMP> inline typename Modular<uint8_t, COMP>::Element&
    Modular<uint8_t, COMP>::init (Element& r, const uint32_t a) const
    {
	r = static_cast<Element>(a % uint32_t(_p));
	return r;
    }

    template<typename COMP> inline typename Modular<uint8_t, COMP>::Element&
    Modular<uint8_t, COMP>::init (Element& r, const int64_t a) const
    {
	r = static_cast<Element>(std::abs(a) % int64_t(_p));
        if (a < 0) negin(r);
	return r;
    }

    template<typename COMP> inline typename Modular<uint8_t, COMP>::Element&
    Modular<uint8_t, COMP>::init (Element& r, const uint64_t a) const
    {
	r = static_cast<Element>(a % uint64_t(_p));
	return r;
    }

    template<typename COMP> inline typename Modular<uint8_t, COMP>::Element&
    Modular<uint8_t, COMP>::init (Element& r, const Integer& a) const
    {
	r = static_cast<Element>(((a < 0)? -a : a) % _p);
	return r = (a < 0)? negin(r) : r;
    }

    // ------------------------
    // ----- Classic arithmetic

    template<typename COMP>
    inline typename Modular<uint8_t, COMP>::Element &Modular<uint8_t, COMP>::add
    (Element &x, const Element &y, const Element &z) const
    {
        __GIVARO_MODULAR_INTEGER_ADD(x,_p,y,z);
        return x;
    }

    template<typename COMP>
    inline typename Modular<uint8_t, COMP>::Element &Modular<uint8_t, COMP>::sub
    (Element &x, const Element &y, const Element &z) const
    {
        return __GIVARO_MODULAR_INTEGER_SUB(x,_p,y,z);
    }

    template<typename COMP>
    inline typename Modular<uint8_t, COMP>::Element &Modular<uint8_t, COMP>::mul
    (Element &x, const Element &y, const Element &z) const
    {
        return __GIVARO_MODULAR_INTEGER_MUL(x,_p,y,z);
    }

    template<typename COMP>
    inline typename Modular<uint8_t, COMP>::Element &Modular<uint8_t, COMP>::div
    (Element &x, const Element &y, const Element &z) const
    {
        return mulin(inv(x, z), y);
    }

    template<typename COMP>
    inline typename Modular<uint8_t, COMP>::Element &Modular<uint8_t, COMP>::neg
    (Element &x, const Element &y) const
    {
        return __GIVARO_MODULAR_INTEGER_NEG(x,_p,y);
    }

    template<typename COMP>
    inline typename Modular<uint8_t, COMP>::Element &Modular<uint8_t, COMP>::inv
    (Element &x, const Element &y) const
    {
        invext(x, y, _p);
        return (int8_t(x) < 0)? x = Element(x + _p) : x;
    }

    template<typename COMP>
    inline typename Modular<uint8_t, COMP>::Element &Modular<uint8_t, COMP>::addin
    (Element &x, const Element &y) const
    {
        __GIVARO_MODULAR_INTEGER_ADDIN(x,_p,y);
        return x;
    }

    template<typename COMP>
    inline typename Modular<uint8_t, COMP>::Element &Modular<uint8_t, COMP>::subin
    (Element &x, const Element &y) const
    {
        __GIVARO_MODULAR_INTEGER_SUBIN(x,_p,y);
        return x;
    }

    template<typename COMP>
    inline typename Modular<uint8_t, COMP>::Element &Modular<uint8_t, COMP>::mulin
    (Element &x, const Element &y) const
    {
        return __GIVARO_MODULAR_INTEGER_MULIN(x,_p,y);
    }

    template<typename COMP>
    inline typename Modular<uint8_t, COMP>::Element &Modular<uint8_t, COMP>::divin
    (Element &x, const Element &y) const
    {
        typename Modular<uint8_t, COMP>::Element iy;
        return mulin(x, inv(iy, y));
    }

    template<typename COMP>
    inline typename Modular<uint8_t, COMP>::Element &Modular<uint8_t, COMP>::negin
    (Element &x) const
    {
        return __GIVARO_MODULAR_INTEGER_NEGIN(x,_p);
    }

    template<typename COMP>
    inline typename Modular<uint8_t, COMP>::Element &Modular<uint8_t, COMP>::invin
    (Element &x) const
    {
        return inv(x, x);
    }

    // -- axpy: r <- a * x + y
    template<typename COMP>
    inline typename Modular<uint8_t, COMP>::Element &Modular<uint8_t, COMP>::axpy
    (Element &r, const Element &a, const Element &x, const Element &y) const
    {
        __GIVARO_MODULAR_INTEGER_MULADD(r, _p, a, x, y);
        return r;
    }

    template<typename COMP>
    inline typename Modular<uint8_t, COMP>::Element &Modular<uint8_t, COMP>::axpyin
    (Element &r, const Element &a, const Element &x) const
    {
        __GIVARO_MODULAR_INTEGER_MULADDIN(r, _p, a, x);
        return r;
    }

    // -- axmy: r <- a * x - y
    template<typename COMP>
    inline typename Modular<uint8_t, COMP>::Element &Modular<uint8_t, COMP>::axmy
    (Element& r, const Element &a, const Element &x, const Element &y) const
    {
        __GIVARO_MODULAR_INTEGER_MULSUB(r, _p, a, x, y);
        return r;
    }

    template<typename COMP>
    inline typename Modular<uint8_t, COMP>::Element &Modular<uint8_t, COMP>::axmyin
    (Element& r, const Element &a, const Element &x) const
    {
        maxpyin(r,a,x);
        return negin(r);
    }

    // -- maxpy:   r <- y - a * x
    template<typename COMP>
    inline typename Modular<uint8_t, COMP>::Element& Modular<uint8_t, COMP>::maxpy
    (Element& r, const Element& a, const Element& x, const Element& y) const
    {
        r = y;
        __GIVARO_MODULAR_INTEGER_SUBMULIN(r, _p, a, x);
        return r;
    }

    template<typename COMP>
    inline typename Modular<uint8_t, COMP>::Element& Modular<uint8_t, COMP>::maxpyin
    (Element& r, const Element& a, const Element& x) const
    {
        __GIVARO_MODULAR_INTEGER_SUBMULIN(r, _p, a, x);
        return r;
    }

    // ----------------
    // ----- IO methods

    template<>
    inline std::ostream& Modular<uint8_t, int8_t>::write (std::ostream& s) const
    {
        return s << "Modular<uint8_t, uint8_t> modulo " << residu();
    }

    template<>
    inline std::ostream& Modular<uint8_t, uint8_t>::write (std::ostream& s) const
    {
        return s << "Modular<uint8_t, uint8_t> modulo " << residu();
    }

    template<>
    inline std::ostream& Modular<uint8_t, int16_t>::write (std::ostream& s) const
    {
        return s << "Modular<uint8_t, uint16_t> modulo " << residu();
    }

    template<>
    inline std::ostream& Modular<uint8_t, uint16_t>::write (std::ostream& s) const
    {
        return s << "Modular<uint8_t, uint16_t> modulo " << residu();
    }

    template<typename COMP>
    inline std::istream &Modular<uint8_t, COMP>::read (std::istream &is)
    {
        is >> _p;
        return is;
    }

    template<typename COMP>
    inline std::ostream &Modular<uint8_t, COMP>::write (std::ostream &os, const Element &x) const
    {
        return os << int32_t(x);
    }

    template<typename COMP>
    inline std::istream &Modular<uint8_t, COMP>::read (std::istream &is, Element &x) const
    {
        int64_t tmp;
        is >> tmp;
        init(x,tmp);
        return is;
    }

} // namespace Givaro

#endif

