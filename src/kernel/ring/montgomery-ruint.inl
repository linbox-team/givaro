// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: A. Breust (adapted)
// ==========================================================================

#ifndef __GIVARO_montgomery_ruint_INL
#define __GIVARO_montgomery_ruint_INL

namespace Givaro
{
    // ------------------------
    // ----- Internal reduction

    template<size_t K>
    inline typename Montgomery<RecInt::ruint<K>>::Element&
    Montgomery<RecInt::ruint<K>>::mg_reduc(Element& a, const Element& b) const
    {
        bool r;
        Element b0;

        RecInt::mul(b0, b, _p1);                // m = b * p1 mod r
        RecInt::laddmul(r, a, b0, b0, _p, b);   // a|b0 = b * (p1 * p + 1)

        if (r || a >= _p) RecInt::sub(a, _p);
        return a;
    }

    template<size_t K>
    inline typename Montgomery<RecInt::ruint<K>>::Element&
    Montgomery<RecInt::ruint<K>>::mg_reduc(Element& a, const LargeElement& b) const
    {
        bool r;
        Element b0;

        RecInt::mul(b0, b.Low, _p1);            // m = b * p1 mod r
        RecInt::laddmul(r, a, b0, b0, _p, b);   // a|b0 = b * (p1 * p + 1)

        if (r || a >= _p) RecInt::sub(a, _p);
        return a;
    }

    template <size_t K>
    inline typename Montgomery<RecInt::ruint<K>>::Element&
    Montgomery<RecInt::ruint<K>>::to_mg(Element& a, const Element& b) const
    {
        mul(a, b, _r2);
        return a;
    }

    template <size_t K>
    inline typename Montgomery<RecInt::ruint<K>>::Element&
    Montgomery<RecInt::ruint<K>>::to_mg(Element& a) const
    {
        return to_mg(a, a);
    }

    // ------------------------
    // ----- Classic arithmetic

    template<size_t K>
    inline typename Montgomery<RecInt::ruint<K>>::Element& Montgomery<RecInt::ruint<K>>::mul
    (Element& r, const Element& a, const Element& b) const
    {
        LargeElement res;
        RecInt::lmul(res, a, b);
        mg_reduc(r, res);
        return r;
    }

    template<size_t K>
    inline typename Montgomery<RecInt::ruint<K>>::Element& Montgomery<RecInt::ruint<K>>::sub
    (Element& r, const Element& a, const Element& b) const
    {
        if (a < b) { // b > 0 and (a - b) < p
            RecInt::sub(r, _p, b);
            RecInt::add(r, a);
        } else {
            RecInt::sub(r, a, b);
        }
        return r;
    }

    template<size_t K>
    inline typename Montgomery<RecInt::ruint<K>>::Element& Montgomery<RecInt::ruint<K>>::add
    (Element& r, const Element& a, const Element& b) const
    {
        bool ret;
        RecInt::add(ret, r, a, b);
        if (ret || r >= _p) RecInt::sub(r, _p);
        return r;
    }

    template<size_t K>
    inline typename Montgomery<RecInt::ruint<K>>::Element& Montgomery<RecInt::ruint<K>>::neg
    (Element& r, const Element& a) const
    {
        if (a == 0) RecInt::reset(r);
        else        RecInt::sub(r, _p, a);
        return r;
    }

    template<size_t K>
    inline typename Montgomery<RecInt::ruint<K>>::Element& Montgomery<RecInt::ruint<K>>::inv
    (Element& r, const Element& a) const
    {
        RecInt::inv_mod(r, a, _p);  // r = (aR)^-1 mod p
        return mulin(r, _r3);       // r = a^-1 * R mod p
    }

    template<size_t K>
    inline bool Montgomery<RecInt::ruint<K>>::isUnit(const Element& a) const
    {
        RecInt::ruint<K> d;
        gcd(d,a,_p);
        return (d==1) || (d==-1);
    }

    template<size_t K>
    inline typename Montgomery<RecInt::ruint<K>>::Element& Montgomery<RecInt::ruint<K>>::div
    (Element& r, const Element& a, const Element& b) const
    {
        return mulin(inv(r, b), a);
    }

    template<size_t K>
    inline typename Montgomery<RecInt::ruint<K>>::Element& Montgomery<RecInt::ruint<K>>::mulin
    (Element& r, const Element& a) const
    {
        LargeElement res;
        RecInt::lmul(res, r, a);
        return mg_reduc(r, res);
    }

    template<size_t K>
    inline typename Montgomery<RecInt::ruint<K>>::Element& Montgomery<RecInt::ruint<K>>::divin
    (Element& r, const Element& a) const
    {
        Element ia;
        return mulin(r, inv(ia, a));
    }

    template<size_t K>
    inline typename Montgomery<RecInt::ruint<K>>::Element& Montgomery<RecInt::ruint<K>>::addin
    (Element& r, const Element& a) const
    {
        bool ret;
        RecInt::add(ret, r, a);
        if (ret || r >= _p) RecInt::sub(r, _p);
        return r;
    }

    template<size_t K>
    inline typename Montgomery<RecInt::ruint<K>>::Element& Montgomery<RecInt::ruint<K>>::subin
    (Element& r, const Element& a) const
    {
        if (r < a) {
            RecInt::add(r, _p - a);
        } else {
            RecInt::sub(r, a);
        }
        return r;
    }

    template<size_t K>
    inline typename Montgomery<RecInt::ruint<K>>::Element& Montgomery<RecInt::ruint<K>>::negin
    (Element& r) const
    {
        if (r == 0) RecInt::reset(r);
        else        RecInt::sub(r, _p, r);
        return r;
    }
    template<size_t K>
    inline typename Montgomery<RecInt::ruint<K>>::Element& Montgomery<RecInt::ruint<K>>::invin
    (Element& r) const
    {
        return inv(r, r);
    }

    template<size_t K>
    inline typename Montgomery<RecInt::ruint<K>>::Element& Montgomery<RecInt::ruint<K>>::axpy
    (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        mul(r, a, b);
        return addin(r, c);
    }

    template<size_t K>
    inline typename Montgomery<RecInt::ruint<K>>::Element&  Montgomery<RecInt::ruint<K>>::axpyin
    (Element& r, const Element& a, const Element& b) const
    {
        Element res;
        mul(res, a, b);
        return addin(r, res);
    }

    template<size_t K>
    inline typename Montgomery<RecInt::ruint<K>>::Element& Montgomery<RecInt::ruint<K>>::maxpy
    (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        mul(r, a, b);
        return sub(r, c, r);
    }

    template<size_t K>
    inline typename Montgomery<RecInt::ruint<K>>::Element&  Montgomery<RecInt::ruint<K>>::maxpyin
    (Element& r, const Element& a, const Element& b) const
    {
        Element res;
        mul(res, a, b);
        return subin(r, res);
    }

    template<size_t K>
    inline typename Montgomery<RecInt::ruint<K>>::Element&  Montgomery<RecInt::ruint<K>>::axmy
    (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        mul(r, a, b);
        subin(r, c);
        return r;
    }

    template<size_t K>
    inline typename Montgomery<RecInt::ruint<K>>::Element&  Montgomery<RecInt::ruint<K>>::axmyin
    (Element& r, const Element& a, const Element& b) const
    {
        Element res;
        mul(res, a, b);
        return sub(r, res, r);
    }

    //----- IO

    template<size_t K>
    inline std::ostream& Montgomery<RecInt::ruint<K>>::write (std::ostream& s) const
    {
        return s << "Montgomery<RecInt::ruint<" << K << ">> modulo " << residu();
    }

    template<size_t K>
    inline std::istream& Montgomery<RecInt::ruint<K>>::read (std::istream& s, Element& a) const
    {
        Integer tmp;
        s >> tmp;
        init(a, tmp);
        return s;
    }

    template<size_t K>
    inline std::ostream& Montgomery<RecInt::ruint<K>>::write (std::ostream& s, const Element& a) const
    {
        Element ar;
        return s << mg_reduc(ar, a);
    }
}

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
