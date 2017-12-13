// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: A. Breust (adapted)
// ==========================================================================

#ifndef __GIVARO_mod_ruint_INL
#define __GIVARO_mod_ruint_INL

#include "modular-defines.h"

namespace Givaro
{
#define TMPL template<typename Storage_t, typename Compute_t>
#define MOD Mod<Storage_t, Compute_t, \
  typename std::enable_if<is_same_RecInt<Storage_t, Compute_t>::value || \
                          is_smaller_RecInt<Storage_t, Compute_t>::value>::type>

// Storage_t and Compute_t of same size:
#define MOD_SAME Mod<Storage_t, Compute_t, \
  typename std::enable_if<is_same_RecInt<Storage_t, Compute_t>::value>::type>

// Storage_t is ruint<K> and Compute_t is ruint<K+1>
#define MOD_DIFF \
  Mod<Storage_t, Compute_t, typename std::enable_if<is_smaller_RecInt<Storage_t, Compute_t>::value>::type>


    // ------------------------
    // ----- Classic arithmetic
    //TMPL
    //inline
    //typename Mod<Storage_t, Compute_t, typename std::enable_if<is_same_RecInt<Storage_t, Compute_t>::value>::type>::Element&
    //Mod<Storage_t, Compute_t, typename std::enable_if<is_same_RecInt<Storage_t, Compute_t>::value>::type>::mul
    //(Element& r, const Element& a, const Element& b) const
    //{
    //    __GIVARO_MODULAR_RECINT_MUL(r,_p,a,b);
    //    return r;
    //}
    //TMPL
    //inline 
    //typename std::enable_if<is_smaller_RecInt<Storage_t, Compute_t>::value, typename MOD::Element&>::type
    //MOD::mul
    //    (Element& r, const Element& a, const Element& b) const
    //{
    //    __GIVARO_MODULAR_RECINT_LMUL(r,_p,a,b);
    //    return r;
    //}

    //TMPL
    ////inline
    ////typename std::enable_if<is_same_RecInt<Storage_t, Compute_t>::value, typename MOD::Element&>::type
    ////typename MOD::Element&
    //Storage_t&
    //MOD::mul
    //(Element& r, const Element& a, const Element& b) const
    //{
    //    __GIVARO_MODULAR_RECINT_MUL(r,_p,a,b);
    //    return r;
    //}

    TMPL
    inline
    typename MOD::Element&
    MOD::mul
    (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_MODULAR_RECINT_MUL(r,_p,a,b);
        return r;
    }

    //TMPL
    //inline 
    //typename std::enable_if<is_smaller_RecInt<Storage_t, Compute_t>::value, typename MOD::Element&>::type
    //MOD::mul
    //    (Element& r, const Element& a, const Element& b) const
    //{
    //    __GIVARO_MODULAR_RECINT_LMUL(r,_p,a,b);
    //    return r;
    //}

    TMPL
    inline typename MOD::Element& MOD::sub
        (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_MODULAR_RECINT_SUB(r,_p,a,b);
        return r;
    }
 
    TMPL
    inline typename MOD::Element& MOD::add
        (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_MODULAR_RECINT_ADD(r,_p,a,b);
        return r;
    }

    TMPL
    inline typename MOD::Element& MOD::neg
        (Element& r, const Element& a) const
    {
        __GIVARO_MODULAR_RECINT_NEG(r,_p,a);
        return r;
    }

    TMPL
    inline typename MOD::Element& MOD::inv
        (Element& r, const Element& a) const
    {
        return inv_mod(r, a, _p);
    }

    //TMPL
    //inline bool MOD::isUnit(const Element& a) const 
    //{ 
    //    Element d; 
    //    gcd(d,a,_p); 
    //    return isOne(d) || isMOne(d); 
    //}

    TMPL
    inline typename MOD::Element& MOD::div
        (Element& r, const Element& a, const Element& b) const
    {
        return mulin( inv(r,b), a );
    }

    TMPL
    //inline 
    typename MOD_SAME::Element& MOD_SAME::mulin
        (Element& r, const Element& a) const
    {
        __GIVARO_MODULAR_RECINT_MULIN(r,_p,a);
        return r;
    }
    TMPL
    inline typename MOD_DIFF::Element& MOD_DIFF::mulin
        (Element& r, const Element& a) const
    {
        __GIVARO_MODULAR_RECINT_LMULIN(r,_p,a);
        return r;
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
        __GIVARO_MODULAR_RECINT_ADDIN(r,_p,a);
        return r;
    }

    TMPL
    inline typename MOD::Element& MOD::subin
        (Element& r, const Element& a) const
    {
        __GIVARO_MODULAR_RECINT_SUBIN(r,_p,a);
        return r;
    }

    TMPL
    inline typename MOD::Element& MOD::negin
        (Element& r) const
    {
        __GIVARO_MODULAR_RECINT_NEGIN(r,_p);
        return r;
    }

    TMPL
    inline typename MOD::Element& MOD::invin
        (Element& r) const
    {
        return inv(r, r);
    }

    //TMPL
    //inline typename MOD_SAME::Element& MOD_SAME::axpy
    //    (Element& r, const Element& a, const Element& b, const Element& c) const
    //{
    //    __GIVARO_MODULAR_RECINT_MULADD(r,_p,a,b,c);
    //    return r;
    //}
    TMPL
    inline typename MOD::Element& MOD::axpy
        (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        __GIVARO_MODULAR_RECINT_LMULADD(r,_p,a,b,c);
        return r;
    }

    //TMPL
    //inline typename MOD_SAME::Element&  MOD_SAME::axpyin
    //    (Element& r, const Element& a, const Element& b) const
    //{
    //    __GIVARO_MODULAR_RECINT_MULADDIN(r,_p,a,b);
    //    return r;
    //}
    TMPL
    inline typename MOD::Element&  MOD::axpyin
        (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_MODULAR_RECINT_LMULADDIN(r,_p,a,b);
        return r;
    }

    //TMPL
    //inline typename MOD_SAME::Element& MOD_SAME::maxpy
    //    (Element& r, const Element& a, const Element& b, const Element& c) const
    //{
    //    __GIVARO_MODULAR_RECINT_MUL(r,_p,a,b);
    //    __GIVARO_MODULAR_RECINT_SUB(r,_p,c,r);
    //    return r;
    //}
    TMPL
    inline typename MOD::Element& MOD::maxpy
        (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        __GIVARO_MODULAR_RECINT_LMUL(r,_p,a,b);
        __GIVARO_MODULAR_RECINT_SUB(r,_p,c,r);
        return r;
    }

    //TMPL
    //inline typename MOD_SAME::Element&  MOD_SAME::axmy
    //    (Element& r, const Element& a, const Element& b, const Element& c) const
    //{
    //    __GIVARO_MODULAR_RECINT_MULSUB(r,_p,a,b,c);
    //    return r;
    //}
    TMPL
    inline typename MOD::Element&  MOD::axmy
        (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        __GIVARO_MODULAR_RECINT_LMULSUB(r,_p,a,b,c);
        return r;
    }

    //TMPL
    //inline typename MOD_SAME::Element&  MOD_SAME::maxpyin
    //    (Element& r, const Element& a, const Element& b) const
    //{
    //    __GIVARO_MODULAR_RECINT_SUBMULIN(r,_p,a,b);
    //    __GIVARO_MODULAR_RECINT_NEGIN(r,_p);
    //    return r;
    //}
    TMPL
    inline typename MOD::Element&  MOD::maxpyin
        (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_MODULAR_RECINT_LSUBMULIN(r,_p,a,b);
        __GIVARO_MODULAR_RECINT_NEGIN(r,_p);
        return r;
    }

    //TMPL
    //inline typename MOD_SAME::Element&  MOD_SAME::axmyin
    //    (Element& r, const Element& a, const Element& b) const
    //{
    //    Element rc(r);
    //    __GIVARO_MODULAR_RECINT_MULSUB(r,_p,a,b,rc);
    //    return r;
    //}
    TMPL
    inline typename MOD::Element&  MOD::axmyin
        (Element& r, const Element& a, const Element& b) const
    {
        Element rc(r);
        __GIVARO_MODULAR_RECINT_LMULSUB(r,_p,a,b,rc);
        return r;
    }

#undef TMPL
#undef MOD
#undef MOD_SAME
#undef MOD_DIFF

}


#endif
