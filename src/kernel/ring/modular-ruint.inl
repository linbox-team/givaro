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
        // ----- Modular<ruint<K> >

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K>>::Residu_t
    Modular<RecInt::ruint<K>, RecInt::ruint<K>>::maxCardinality() {
        Residu_t max; max.High=1U; return max;
    }

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
    inline bool Modular<RecInt::ruint<K>, RecInt::ruint<K>>::isUnit(const Element& a) const 
    { 
        Element d; 
        gcd(d,a,_p); 
        return isOne(d) || isMOne(d); 
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
    inline std::ostream& Modular<RecInt::ruint<K>, RecInt::ruint<K>>::write (std::ostream& s, const Element& a) const
    {
        return s << a;
    }  
}

    // -------------
    // ----- Modular<ruint<K>, ruint<K+1> >
namespace Givaro {

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::Residu_t
    Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::maxCardinality() {
        Residu_t max; return RecInt::fill_with_1(max)>>1;
    }

    // ------------------------
    // ----- Classic arithmetic

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::mul
        (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_MODULAR_RECINT_LMUL(r,_p,a,b);
        return r;
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::sub
        (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_MODULAR_RECINT_SUB(r,_p,a,b);
        return r;
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::add
        (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_MODULAR_RECINT_ADD(r,_p,a,b);
        return r;
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::neg
        (Element& r, const Element& a) const
    {
        __GIVARO_MODULAR_RECINT_NEG(r,_p,a);
        return r;
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::inv
        (Element& r, const Element& a) const
    {
        return inv_mod(r, a, _p);
    }

    template<size_t K>
    inline bool Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::isUnit(const Element& a) const 
    { 
        Element d; 
        gcd(d,a,_p); 
        return isOne(d) || isMOne(d); 
    }

  
    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::div
        (Element& r, const Element& a, const Element& b) const
    {
        return mulin( inv(r,b), a );
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::mulin
        (Element& r, const Element& a) const
    {
        __GIVARO_MODULAR_RECINT_LMULIN(r,_p,a);
        return r;
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::divin
        (Element& r, const Element& a) const
    {
        Element ia;
        return mulin(r, inv(ia, a));
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::addin
        (Element& r, const Element& a) const
    {
        __GIVARO_MODULAR_RECINT_ADDIN(r,_p,a);
        return r;
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::subin
        (Element& r, const Element& a) const
    {
        __GIVARO_MODULAR_RECINT_SUBIN(r,_p,a);
        return r;
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::negin
        (Element& r) const
    {
        __GIVARO_MODULAR_RECINT_NEGIN(r,_p);
        return r;
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::invin
        (Element& r) const
    {
        return inv(r, r);
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::axpy
        (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        __GIVARO_MODULAR_RECINT_LMULADD(r,_p,a,b,c);
        return r;
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::Element&  Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::axpyin
        (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_MODULAR_RECINT_LMULADDIN(r,_p,a,b);
        return r;
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::Element& Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::maxpy
        (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        __GIVARO_MODULAR_RECINT_LMUL(r,_p,a,b);
        __GIVARO_MODULAR_RECINT_SUB(r,_p,c,r);
        return r;
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::Element&  Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::axmy
        (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        __GIVARO_MODULAR_RECINT_LMULSUB(r,_p,a,b,c);
        return r;
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::Element&  Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::maxpyin
        (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_MODULAR_RECINT_LSUBMULIN(r,_p,a,b);
        __GIVARO_MODULAR_RECINT_NEGIN(r,_p);
        return r;
    }

    template<size_t K>
    inline typename Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::Element&  Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::axmyin
        (Element& r, const Element& a, const Element& b) const
    {
        Element rc(r);
        __GIVARO_MODULAR_RECINT_LMULSUB(r,_p,a,b,rc);
        return r;
    }

    //----- IO
    
    template<size_t K>
    inline std::ostream& Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::write (std::ostream& s) const
    {
        return s << "Modular<RecInt::ruint<" << K << ">, RecInt::ruint<" << K+1 << ">> modulo " << residu();
    }

    template<size_t K>
    inline std::istream& Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::read (std::istream& s, typename Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::Element& a) const
    {
        Integer tmp;
        s >> tmp;
        init(a, tmp);
        return s;
    }

    template<size_t K>
    inline std::ostream& Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::write (std::ostream& s, const typename Modular<RecInt::ruint<K>, RecInt::ruint<K+1>>::Element& a) const
    {
        return s << a;
    }

} // end of namespace Givaro

#endif
