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

#include "modular-defines.h"

namespace Givaro
{
    // -------------
    // ----- Modular

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Residu_t
        Modular<RecInt::ruint<K>, RecInt::ruint<K>>::maxCardinality() { return RecInt::ruint<K>::maxCardinality(); }

    // ------------------------
    // ----- Classic arithmetic

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::mul
        (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_MODULAR_RECINT_MUL(r,_p,a,b);
        return r;
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::sub
        (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_MODULAR_RECINT_SUB(r,_p,a,b);
        return r;
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::add
        (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_MODULAR_RECINT_ADD(r,_p,a,b);
        return r;
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::neg
        (Element& r, const Element& a) const
    {
        __GIVARO_MODULAR_RECINT_NEG(r,_p,a);
        return r;
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::inv
        (Element& r, const Element& a) const
    {
        return inv_mod(r, a, _p);
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::div
        (Element& r, const Element& a, const Element& b) const
    {
        return mulin( inv(r,b), a );
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::mulin
        (Element& r, const Element& a) const
    {
        __GIVARO_MODULAR_RECINT_MULIN(r,_p,a);
        return r;
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::divin
        (Element& r, const Element& a) const
    {
        Element ia;
        return mulin(r, inv(ia, a));
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::addin
        (Element& r, const Element& a) const
    {
        __GIVARO_MODULAR_RECINT_ADDIN(r,_p,a);
        return r;
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::subin
        (Element& r, const Element& a) const
    {
        __GIVARO_MODULAR_RECINT_SUBIN(r,_p,a);
        return r;
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::negin
        (Element& r) const
    {
        __GIVARO_MODULAR_RECINT_NEGIN(r,_p);
        return r;
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::invin
        (Element& r) const
    {
        return inv(r, r);
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::axpy
        (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        __GIVARO_MODULAR_RECINT_MULADD(r,_p,a,b,c);
        return r;
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element&  Modular<RecInt::ruint<K>, RecInt::ruint<K>>::axpyin
        (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_MODULAR_RECINT_MULADDIN(r,_p,a,b);
        return r;
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::maxpy
        (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        __GIVARO_MODULAR_RECINT_MUL(r,_p,a,b);
        __GIVARO_MODULAR_RECINT_SUB(r,_p,c,r);
        return r;
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element&  Modular<RecInt::ruint<K>, RecInt::ruint<K>>::axmy
        (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        __GIVARO_MODULAR_RECINT_MULSUB(r,_p,a,b,c);
        return r;
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element&  Modular<RecInt::ruint<K>, RecInt::ruint<K>>::maxpyin
        (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_MODULAR_RECINT_SUBMULIN(r,_p,a,b);
        __GIVARO_MODULAR_RECINT_NEGIN(r,_p);
        return r;
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Element&  Modular<RecInt::ruint<K>, RecInt::ruint<K>>::axmyin
        (Element& r, const Element& a, const Element& b) const
    {
        Element rc(r);
        __GIVARO_MODULAR_RECINT_MULSUB(r,_p,a,b,rc);
        return r;
    }

    //----- IO
    
    template<size_t K>
    inline std::ostream& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::write (std::ostream& s) const
    {
        return s << "Modular<RecInt::ruint<" << K << ">, RecInt::ruint<" << K << ">> modulo " << residu();
    }

    template<size_t K>
    inline std::istream& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::read (std::istream& s, Element& a) const
    {
        Integer tmp;
        s >> tmp;
        init(a, tmp);
        return s;
    }

    template<size_t K>
    inline std::ostream& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::write (std::ostream& s, const Element a) const
    {
        return s << a;
    }
}

#endif
