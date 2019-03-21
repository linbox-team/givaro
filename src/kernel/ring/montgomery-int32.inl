// ==========================================================================
// $Id: givmontg32.inl,v 1.14 2011-02-04 14:11:46 jgdumas Exp $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// ==========================================================================

#ifndef __GIVARO_mong32_INL
#define __GIVARO_mong32_INL

namespace Givaro {

    inline Montgomery<int32_t>::Element Montgomery<int32_t>::redcal(const Element c) const
    {
        Element c0 = (Element)(c & MASK32);	/* c mod B */
        c0 = (Element)((c0 * _nim) & MASK32); 	/* -c/p mod B */
        // c0 *= _p;
        // c0 += c;			/* c = 0 mod B */
        c0 = c + c0*_p;
        c0 >>= HALF_BITS32;
        return (c0>=_p?c0-=_p:c0);
    }
    inline Montgomery<int32_t>::Element Montgomery<int32_t>::redcsal(const Element c) const
    {
        Element c0 = (Element)((c * _nim) & MASK32); 	/* -c/p mod B */
        c0 = c + c0 * _p; 		/* c = 0 mod B */
        c0 >>= HALF_BITS32;
        return (c0>=_p?c0-=_p:c0);
    }

    inline Montgomery<int32_t>::Element& Montgomery<int32_t>::redc(Element& r, const Element c) const
    {
        r = (Element)(c & MASK32);			/* c mod B */
        r *= _nim;
        r &= MASK32;
        r *= _p;
        r += c;
        r >>= HALF_BITS32;
        return (r>=_p?r-=_p:r);
    }

    inline Montgomery<int32_t>::Element& Montgomery<int32_t>::redcs(Element& r, const Element c) const
    {
        r = (Element)((c * _nim) & MASK32); 	/* -c/p mod B */
        r = c + r * _p; 		/* c = 0 mod B */
        r >>= HALF_BITS32;
        return (r>=_p?r-=_p:r);
    }

    inline Montgomery<int32_t>::Element& Montgomery<int32_t>::redcin(Element& r) const
    {
        Element c0 = (Element)(r & MASK32);	/* c mod B */
        c0 = (Element)((c0 * _nim) & MASK32); 	/* -c/p mod B */
        r += c0 * _p; 			/* c = 0 mod B */
        r >>= HALF_BITS32;
        return (r>=_p?r-=_p:r);
    }
    inline Montgomery<int32_t>::Element& Montgomery<int32_t>::redcsin(Element& r) const
    {
        Element c0 = (Element)((r * _nim) & MASK32); 	/* -c/p mod B */
        r += c0 * _p; 				/* c = 0 mod B */
        r >>= HALF_BITS32;
        return (r>=_p?r-=_p:r);
    }

} // namespace Givaro

// r = a*b
#define __GIVARO_MONTG32_MUL(r,p,a,b) (redc(r,a*b))
// r *= a
#define __GIVARO_MONTG32_MULIN(r,p,a) (redcin(r*=a))

// r = a - b
#define __GIVARO_MONTG32_SUB(r,p,a,b) ( r = (a>=b)? a-b: (p-b)+a )
// r -= a
//#define __GIVARO_MONTG32_SUBIN(r,p,a) { r -= a; r= (r < 0 ? r+p : r); }
#define __GIVARO_MONTG32_SUBIN(r,p,a) { if (r<a) r+=(p-a); else r-=a; }

// r = a+b
#define __GIVARO_MONTG32_ADD(r,p,a,b) { r = (a+b); r= (r < p ? r : r-p); }
// r += a
#define __GIVARO_MONTG32_ADDIN(r,p,a) { r += a;  r= (r < p ? r : r-p); }

// r <- a*b+c % p
//#define __GIVARO_MONTG32_MULADD(r,p,a,b,c) { redc(r,a*b) += c; r= (r < p ? r : r-p); }
#define __GIVARO_MONTG32_MULADD(r,p,a,b,c) { r = redcal(a*b) + c; r= (r < p ? r : r-p); }


#define __GIVARO_MONTG32_MULADDIN(r,p,a,b) { r+=redcal(a*b); r= (r < p ? r : r-p); }

#define __GIVARO_MONTG32_NEG(r,p,a) (r = (a == 0 ? 0 : p-a))
#define __GIVARO_MONTG32_NEGIN(r,p) (r = (r == 0 ? 0 : p-r))

namespace Givaro {

    inline Montgomery<int32_t>::Element&
    Montgomery<int32_t>::mul (Element& r, const Element& a, const Element& b) const
    {
        return __GIVARO_MONTG32_MUL(r,_p,a,b);
    }

    inline Montgomery<int32_t>::Element&
    Montgomery<int32_t>::sub (Element& r, const Element& a, const Element& b) const
    {
        return __GIVARO_MONTG32_SUB(r,_p,a,b);
    }

    inline Montgomery<int32_t>::Element&
    Montgomery<int32_t>::add (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_MONTG32_ADD(r,_p,a,b);
        return r;
    }

    inline Montgomery<int32_t>::Element&
    Montgomery<int32_t>::neg (Element& r, const Element& a) const
    {
        return __GIVARO_MONTG32_NEG(r,_p,a);
    }

    inline Montgomery<int32_t>::Element&
    Montgomery<int32_t>::inv (Element& r, const Element& a) const
    {
        // invext(aB) --> 1/a*1/B
        // % * B^3    --> B²/a
        // redc       --> B/a
        int32_t t;
        invext(t, int32_t(a), int32_t(_p));
        if (t < 0) t += _p;
        return redc(r, uint32_t(t) * _B3p) ;
    }

    inline bool Montgomery<int32_t>::isUnit(const Element& a) const
    {
        int32_t u,d;
        extended_euclid(u,d,int32_t(a),int32_t(_p));
        return (d==1)||(d==-1);
    }

    inline Montgomery<int32_t>::Element&
    Montgomery<int32_t>::div (Element& r, const Element& a, const Element& b) const
    {
        return mulin( inv(r,b), a );
    }

    inline Montgomery<int32_t>::Element&
    Montgomery<int32_t>::mulin (Element& r, const Element& a) const
    {
        return __GIVARO_MONTG32_MULIN(r,_p, a);
    }

    inline Montgomery<int32_t>::Element&
    Montgomery<int32_t>::divin (Element& r, const Element& a) const
    {
        Montgomery<int32_t>::Element ia;
        inv(ia, a);
        return mulin(r, ia);
    }

    inline Montgomery<int32_t>::Element&
    Montgomery<int32_t>::addin (Element& r, const Element& a) const
    {
        __GIVARO_MONTG32_ADDIN(r,_p, a);
        return r;
    }

    inline Montgomery<int32_t>::Element&
    Montgomery<int32_t>::subin (Element& r, const Element& a) const
    {
        __GIVARO_MONTG32_SUBIN(r,_p, a);
        return r;
    }


    inline Montgomery<int32_t>::Element&
    Montgomery<int32_t>::negin (Element& r) const
    {
        return __GIVARO_MONTG32_NEGIN(r,_p);
    }

    inline Montgomery<int32_t>::Element&
    Montgomery<int32_t>::invin (Element& r) const
    {
        Element t;
        return r = inv(t,r);
    }

    inline Montgomery<int32_t>::Element&
    Montgomery<int32_t>::axpy(Element& r, const Element& a, const Element& b, const Element& c) const
    {
        __GIVARO_MONTG32_MULADD(r, _p, a, b, c);
        return r;
    }

    inline Montgomery<int32_t>::Element&
    Montgomery<int32_t>::axpyin(Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_MONTG32_MULADDIN(r, _p, a, b);
        return r;
    }

    // r <- a*b-c
    inline Montgomery<int32_t>::Element&
    Montgomery<int32_t>::axmy(Element& r, const Element& a, const Element& b, const Element& c) const
    {
        return this->subin(this->mul(r,a,b), c);
    }

    // r = c - a*b
    inline Montgomery<int32_t>::Element&
    Montgomery<int32_t>::maxpy(Element& r, const Element& a, const Element& b, const Element& c) const
    {
        Element t;
        return this->sub(r, c, this->mul(t,a,b));
    }

    // r -= a*b
    inline Montgomery<int32_t>::Element&
    Montgomery<int32_t>::maxpyin(Element& r, const Element& a, const Element& b) const
    {
        Element t;
        return this->subin(r, this->mul(t,a,b));
    }

    // r = a*b - r
    inline Montgomery<int32_t>::Element&
    Montgomery<int32_t>::axmyin(Element& r, const Element& a, const Element& b) const
    {
        maxpyin(r,a,b);
        return negin(r);
    }

    // --------------------
    // ----- Initialisation

    inline Montgomery<int32_t>::Element&
    Montgomery<int32_t>::init (Element& r, const double a) const
    {
        r = static_cast<Element>(std::fmod((a < 0.0)? -a : a, double(_p)));
        if (a < 0.0) negin(r);
        return redc(r, r * _B2p);
    }

    inline Montgomery<int32_t>::Element&
    Montgomery<int32_t>::init (Element& r, const int64_t a) const
    {

        r = static_cast<Element>(std::abs(a) % int64_t(_p));
        if (a < 0) negin(r);
        return redc(r, r * _B2p);
    }

    inline Montgomery<int32_t>::Element&
    Montgomery<int32_t>::init (Element& r, const uint64_t a) const
    {
        r = static_cast<Element>(a % uint64_t(_p));
        return redc(r, r * _B2p);
    }

    inline Montgomery<int32_t>::Element&
    Montgomery<int32_t>::init (Element& r, const Integer& a) const
    {
        r = static_cast<Element>(((a < 0)? -a : a) % _p);
        if (a < 0) negin(r);
        return redc(r, r * _B2p);
    }

    //----- IO

    inline std::ostream&
    Montgomery<int32_t>::write (std::ostream& s) const
    {
        return s << "Montgomery<int32_t> modulo " << residu();
    }

    inline std::istream&
    Montgomery<int32_t>::read (std::istream& s, Element& a) const
    {
        Integer tmp;
        s >> tmp;
        init(a, tmp);
        return s;
    }

    inline std::ostream&
    Montgomery<int32_t>::write (std::ostream& s, const Element& a) const
    {
        //         Element tmp;
        //         return s << redcs(tmp,a);
        Element tmp;
        redcs(tmp,a);
        return s << tmp;
    }

} // namespace Givaro

#endif // __GIVARO_mong32_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
