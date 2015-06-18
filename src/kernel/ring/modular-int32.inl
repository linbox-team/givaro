// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givzpz32std.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givzpz32std.inl,v 1.20 2011-02-02 17:16:43 bboyer Exp $
// ==========================================================================

#ifndef __GIVARO_modular_int32_INL
#define __GIVARO_modular_int32_INL

#include <cmath>

#include "modular-defines.h"

namespace Givaro {

    // -------------
    // ----- Modular 

    template<>
    inline Modular<int32_t, int32_t>::Residu_t
    Modular<int32_t, int32_t>::maxCardinality() { return 65535u; } // 2^16 - 1

    template<>
    inline Modular<int32_t, uint32_t>::Residu_t
    Modular<int32_t, uint32_t>::maxCardinality() { return 65535u; }

    template<>
    inline Modular<int32_t, uint64_t>::Residu_t
    Modular<int32_t, uint64_t>::maxCardinality() { return 2147483647u; } // 2^31 - 1

    template<>
    inline Modular<int32_t, int64_t>::Residu_t
    Modular<int32_t, int64_t>::maxCardinality() { return 2147483647u; }

    // ------------------------
    // ----- Classic arithmetic
        
    template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::mul
    (Element& r, const Element& a, const Element& b) const
    {
        return  __GIVARO_MODULAR_INTEGER_MUL(r,_p,a,b);
    }

    template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::sub
    (Element& r, const Element& a, const Element& b) const
    {
        return __GIVARO_MODULAR_INTEGER_SUB(r,_p,a,b);
    }

    template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::add
    (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_MODULAR_INTEGER_ADD(r,_p,a,b);
        return r;
    }

    template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::neg
    (Element& r, const Element& a) const
    {
        return __GIVARO_MODULAR_INTEGER_NEG(r,_p,a);
    }

    template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::inv
    (Element& r, const Element& a) const
    {
        invext(r, a, int32_t(_p));
        return (r < 0)? r += int32_t(_p) : r;
    }

    template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::div
    (Element& r, const Element& a, const Element& b) const
    {
        return mulin(inv(r, b), a);
    }
    template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::mulin
    (Element& r, const Element& a) const
    {
        r = Element(Compute_t(r)*Compute_t(a) % Compute_t(_p));
        return r;
    }

    template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::divin
    (Element& r, const Element& a) const
    {
        Element ia;
        return mulin(r, inv(ia, a));
    }

    template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::addin
    (Element& r, const Element& a) const
    {
        int32_t tmp = r;
        __GIVARO_MODULAR_INTEGER_ADDIN(tmp,_p, a);
        return r = Element(tmp);
    }

    template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::subin
    (Element& r, const Element& a) const
    {
        int32_t tmp = r;
        __GIVARO_MODULAR_INTEGER_SUBIN(tmp,_p, a);
        return r = (Modular<int32_t, COMP>::Element)tmp;
    }

    template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::negin
    (Element& r) const
    {
        return __GIVARO_MODULAR_INTEGER_NEGIN(r,_p);
    }

    template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::invin
    (Element& r) const
    {
        return inv(r, r);
    }
        
    template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::axpy
    (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        return __GIVARO_MODULAR_INTEGER_MULADD(r, _p, a, b, c);
    }

    template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element&  Modular<int32_t, COMP>::axpyin
    (Element& r, const Element& a, const Element& b) const
    {
        return __GIVARO_MODULAR_INTEGER_MULADDIN(r, _p, a, b);
    }
        
    template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::maxpy
    (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        int32_t tmp;
        __GIVARO_MODULAR_INTEGER_MUL(tmp, _p, a, b);
        __GIVARO_MODULAR_INTEGER_SUB(r, _p, c, tmp);
        return r;
    }

    template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element&  Modular<int32_t, COMP>::axmy
    (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        int32_t tmp;
        __GIVARO_MODULAR_INTEGER_MULSUB(tmp, _p, a, b, c);
        return r = (Modular<int32_t, COMP>::Element)tmp;
    }

    template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element&  Modular<int32_t, COMP>::maxpyin
    (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_MODULAR_INTEGER_SUBMULIN(r, _p, a, b );
        return r;
    }

    template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element&  Modular<int32_t, COMP>::axmyin
    (Element& r, const Element& a, const Element& b) const
    {
        return __GIVARO_MODULAR_INTEGER_MULSUB(r, _p, a, b, r );
    }

    // --------------------
    // ----- Initialisation
    
    template<typename COMP> inline typename Modular<int32_t, COMP>::Element&
    Modular<int32_t, COMP>::init(Element& x) const
    {
        return x = zero;
    }
    
    template<typename COMP> inline typename Modular<int32_t, COMP>::Element&
    Modular<int32_t, COMP>::init (Element& x, const float y) const
    {
        Caster<Element,float>(x, std::fmod(y, float(_p)));
        return (x < 0.f ? Caster(x, (x+_p) ) : x);
//         if (x < 0.f) x = static_cast<Element>(x + _p);
//         return x;
    }
    
    template<typename COMP> inline typename Modular<int32_t, COMP>::Element&
    Modular<int32_t, COMP>::init (Element& x, const double y) const
    {
        Caster<Element,double>(x,std::fmod(y, double(_p)));
        return (x < 0.0 ? Caster(x, (x+_p) ) : x);
//         if (x < 0.0) x = static_cast<Element>(x + _p);
//         return x;
    }
    
    template<typename COMP> inline typename Modular<int32_t, COMP>::Element&
    Modular<int32_t, COMP>::init (Element& x, const int64_t y) const
    {
        Caster<Element,int64_t>(x,y % int64_t(_p));
        return (x < 0 ? Caster(x, (x+_p) ) : x);
//         if (x < 0) x = static_cast<Element>(x + _p);
//         return x;
    }
    
    template<typename COMP> inline typename Modular<int32_t, COMP>::Element&
    Modular<int32_t, COMP>::init (Element& x, const uint64_t y) const
    {
        return Caster<Element,uint64_t>(x,y % uint64_t(_p));
    }
    
    template<typename COMP> inline typename Modular<int32_t, COMP>::Element&
    Modular<int32_t, COMP>::init (Element& x, const Integer& y) const
    {
        Caster<Element,Residu_t>(x, y % _p);
        return (x < 0 ? Caster(x, (x+_p) ) : x);
//         if (x < 0) x = static_cast<Element>(x + _p);
//         return x;
    }

    template<typename COMP> inline  typename Modular<int32_t, COMP>::Element&
    Modular<int32_t, COMP>::assign ( Element& r, const Element& a ) const
    {
        return r = a;
    }
        
    // ------------
    // ----- Reduce

    template<typename COMP> inline typename Modular<int32_t, COMP>::Element&
    Modular<int32_t, COMP>::reduce(Element& r, const Element& a) const
    {
        r = a % static_cast<Element>(_p);
        return (r < 0 ? r = static_cast<Element>(r + _p) : r);
    }

    template<typename COMP> inline typename Modular<int32_t, COMP>::Element&
    Modular<int32_t, COMP>::reduce(Element& r) const
    {
        r %= static_cast<Element>(_p);
        return (r < 0 ? r = static_cast<Element>(r + _p) : r);
    }
        
    // --------
    // ----- IO

    template<>
    inline std::ostream& Modular<int32_t, int32_t>::write (std::ostream& s) const
    {
        return s << "Modular<int32_t, uint32_t> modulo " << residu();
    }

    template<>
    inline std::ostream& Modular<int32_t, uint32_t>::write (std::ostream& s) const
    {
        return s << "Modular<int32_t, uint32_t> modulo " << residu();
    }

    template<>
    inline std::ostream& Modular<int32_t, int64_t>::write (std::ostream& s) const
    {
        return s << "Modular<int32_t, uint64_t> modulo " << residu();
    }

    template<>
    inline std::ostream& Modular<int32_t, uint64_t>::write (std::ostream& s) const
    {
        return s << "Modular<int32_t, uint64_t> modulo " << residu();
    }

    template<typename COMP>
    inline std::istream& Modular<int32_t, COMP>::read (std::istream& s, Element& a) const
    {
        Integer tmp;
        s >> tmp;
        init(a, tmp);
        return s;
    }

    template<typename COMP>
    inline std::ostream& Modular<int32_t, COMP>::write (std::ostream& s, const Element a) const
    {
        return s << a;
    }

} // namespace Givaro

#endif // __GIVARO_zpz32std_INL

