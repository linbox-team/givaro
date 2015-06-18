// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givzpz32uns.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givzpz32uns.inl,v 1.15 2011-02-02 17:16:43 bboyer Exp $
// ==========================================================================

#ifndef __GIVARO_modular_uint32_INL
#define __GIVARO_modular_uint32_INL

#include "givaro/modular-defines.h"
#include <cmath>

namespace Givaro
{

    // -------------
    // ----- Modular 

    template<>
    inline Modular<uint32_t, int32_t>::Residu_t
    Modular<uint32_t, int32_t>::maxCardinality() { return 65535u; } // 2^16 - 1

    template<>
    inline Modular<uint32_t, uint32_t>::Residu_t
    Modular<uint32_t, uint32_t>::maxCardinality() { return 65535u; }

    template<>
    inline Modular<uint32_t, uint64_t>::Residu_t
    Modular<uint32_t, uint64_t>::maxCardinality() { return 2147483647u; } // 2^31 - 1

    template<>
    inline Modular<uint32_t, int64_t>::Residu_t
    Modular<uint32_t, int64_t>::maxCardinality() { return 2147483647u; }

    // ------------------------
    // ----- Classic arithmetic

    template<typename COMP>
    inline typename Modular<uint32_t, COMP>::Element& Modular<uint32_t, COMP>::mul
    (Element& r, const Element& a, const Element& b) const
    {
        return __GIVARO_MODULAR_INTEGER_MUL(r,_p,a,b);
    }

    template<typename COMP>
    inline typename Modular<uint32_t, COMP>::Element& Modular<uint32_t, COMP>::sub
    (Element& r, const Element& a, const Element& b) const
    {
        return __GIVARO_MODULAR_INTEGER_SUB(r,_p,a,b);
    }

    template<typename COMP>
    inline typename Modular<uint32_t, COMP>::Element& Modular<uint32_t, COMP>::add
    (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_MODULAR_INTEGER_ADD(r,_p,a,b);
        return r;
    }

    template<typename COMP>
    inline typename Modular<uint32_t, COMP>::Element& Modular<uint32_t, COMP>::neg
    (Element& r, const Element& a) const
    {
        return __GIVARO_MODULAR_INTEGER_NEG(r,_p,a);
    }

    template<typename COMP>
    inline typename Modular<uint32_t, COMP>::Element& Modular<uint32_t, COMP>::inv
    (Element& r, const Element& a) const
    {
        invext(r, a, _p);
        return (int32_t(r) < 0)? r += _p : r;
    }

    template<typename COMP>
    inline typename Modular<uint32_t, COMP>::Element& Modular<uint32_t, COMP>::div
    (Element& r, const Element& a, const Element& b) const
    {
        return mulin( inv(r, b), a );
    }


    template<typename COMP>
    inline typename Modular<uint32_t, COMP>::Element& Modular<uint32_t, COMP>::mulin
    (Element& r, const Element& a) const
    {
        return __GIVARO_MODULAR_INTEGER_MULIN(r,_p, a);
    }

    template<typename COMP>
    inline typename Modular<uint32_t, COMP>::Element& Modular<uint32_t, COMP>::divin
    (Element& r, const Element& a) const
    {
        typename Modular<uint32_t, COMP>::Element ia;
        inv(ia, a);
        return mulin(r, ia);
    }

    template<typename COMP>
    inline typename Modular<uint32_t, COMP>::Element& Modular<uint32_t, COMP>::addin
    (Element& r, const Element& a) const
    {
        uint32_t tmp = r;
        __GIVARO_MODULAR_INTEGER_ADDIN(tmp,_p, a);
        return r = (Element)tmp;
    }

    template<typename COMP>
    inline typename Modular<uint32_t, COMP>::Element& Modular<uint32_t, COMP>::subin
    (Element& r, const Element& a) const
    {
        uint32_t tmp = r;
        __GIVARO_MODULAR_INTEGER_SUBIN(tmp,_p, a);
        return r = (Element)tmp;
    }


    template<typename COMP>
    inline typename Modular<uint32_t, COMP>::Element& Modular<uint32_t, COMP>::negin
    (Element& r) const
    {
        return __GIVARO_MODULAR_INTEGER_NEGIN(r,_p);
    }

    template<typename COMP>
    inline typename Modular<uint32_t, COMP>::Element& Modular<uint32_t, COMP>::invin
    (Element& r) const
    {
        return inv(r, r);
    }


    template<typename COMP>
    inline typename Modular<uint32_t, COMP>::Element& Modular<uint32_t, COMP>::axpy
    (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        uint32_t tmp;
        __GIVARO_MODULAR_INTEGER_MULADD(tmp, _p, a, b, c);
        return r = (Element)tmp;
    }

    template<typename COMP>
    inline typename Modular<uint32_t, COMP>::Element&  Modular<uint32_t, COMP>::axpyin
    (Element& r, const Element& a, const Element& b) const
    {
        uint32_t tmp = r;
        __GIVARO_MODULAR_INTEGER_MULADDIN(tmp, _p, a, b);
        return r = (Element)tmp;
    }

    template<typename COMP>
    inline typename Modular<uint32_t, COMP>::Element&  Modular<uint32_t, COMP>::axmy
    (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        uint32_t tmp;
        __GIVARO_MODULAR_INTEGER_MULSUB(tmp, _p, a, b, c);
        return r = (Element)tmp;
    }

    // r = c-a*b
    template<typename COMP>
    inline typename Modular<uint32_t, COMP>::Element&  Modular<uint32_t, COMP>::maxpy
    (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        uint32_t tmp=c;
        __GIVARO_MODULAR_INTEGER_SUBMULIN(tmp, _p, a, b);
        return r = (Element)tmp;
    }

    // r -= a*b
    template<typename COMP>
    inline typename Modular<uint32_t, COMP>::Element&  Modular<uint32_t, COMP>::maxpyin
    (Element& r, const Element& a, const Element& b) const
    {
        uint32_t tmp = r;
        __GIVARO_MODULAR_INTEGER_SUBMULIN(tmp, _p, a, b );
        return r = (Element)tmp;
    }

    // r = a*b - r
    template<typename COMP>
    inline typename Modular<uint32_t, COMP>::Element&  Modular<uint32_t, COMP>::axmyin
    (Element& r, const Element& a, const Element& b) const
    {
        maxpyin(r,a,b);
        return negin(r);
    }

    // --------------------
    // ----- Initialisation

    template<typename COMP> inline typename Modular<uint32_t, COMP>::Element&
    Modular<uint32_t, COMP>::init (Element& r, const double a) const
    {
        Caster<Element,double>(r, std::fmod((a < 0.0)? -a : a, double(_p)));
        return (a < 0.0 ? negin(r) : r);
    }

    template<typename COMP> inline typename Modular<uint32_t, COMP>::Element&
    Modular<uint32_t, COMP>::init (Element& r, const int64_t a) const
    {
        Caster<Element,int64_t>(r,std::abs(a) % int64_t(_p));
        return (a < 0 ? negin(r) : r);
    }

    template<typename COMP> inline typename Modular<uint32_t, COMP>::Element&
    Modular<uint32_t, COMP>::init (Element& r, const uint64_t a) const
    {
        return Caster<Element,uint64_t>(r, a % uint64_t(_p));
    }

    template<typename COMP> inline typename Modular<uint32_t, COMP>::Element&
    Modular<uint32_t, COMP>::init (Element& r, const Integer& a) const
    {
        Caster<Element,Integer>(r,((a < 0)? -a : a) % _p);
        return (a < 0 ? negin(r) : r);
    }

    //----- IO

    template<>
    inline std::ostream& Modular<uint32_t, int32_t>::write (std::ostream& s) const
    {
        return s << "Modular<uint32_t, uint32_t> modulo " << residu();
    }

    template<>
    inline std::ostream& Modular<uint32_t, uint32_t>::write (std::ostream& s) const
    {
        return s << "Modular<uint32_t, uint32_t> modulo " << residu();
    }

    template<>
    inline std::ostream& Modular<uint32_t, int64_t>::write (std::ostream& s) const
    {
        return s << "Modular<uint32_t, uint64_t> modulo " << residu();
    }

    template<>
    inline std::ostream& Modular<uint32_t, uint64_t>::write (std::ostream& s) const
    {
        return s << "Modular<uint32_t, uint64_t> modulo " << residu();
    }

    template<typename COMP>
    inline std::istream& Modular<uint32_t, COMP>::read (std::istream& s, Element& a) const
    {
        Integer tmp;
        s >> tmp;
        init(a, tmp);
        return s;
    }

    template<typename COMP>
    inline std::ostream& Modular<uint32_t, COMP>::write (std::ostream& s, const Element a) const
    {
        return s << a;
    }
}

#endif

