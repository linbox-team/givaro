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

// ---------
// -- normalized operations
// ---------

// r = a*b
// #define __GIVARO_ZPZIntType_N_MUL(r,p,a,b) { r = a*b % p; }
#define __GIVARO_ZPZIntType_N_MUL(r,p,a,b) { r = a; r*=b; r %= p; }
// r *= a
//#define __GIVARO_ZPZIntType_N_MULIN(r,p,a) {  r = (r*a % p);  }
#define __GIVARO_ZPZIntType_N_MULIN(r,p,a) {  r *= a; r %= p;  }

// r = a - b
//#define __GIVARO_ZPZIntType_N_SUB(r,p,a,b) { r = (a-b); r= (r < 0 ? r+p : r); }
#define __GIVARO_ZPZIntType_N_SUB(r,p,a,b) { r = (a>=b) ? a-b: (p-b)+a ; }
// r -= a
// #define __GIVARO_ZPZIntType_N_SUBIN(r,p,a) { r -= a; r= (r < 0 ? r+p : r); }
#define __GIVARO_ZPZIntType_N_SUBIN(r,p,a) { if (r<a) r+=(p-a); else r-=a; }

// r = a+b
// #define __GIVARO_ZPZIntType_N_ADD(r,p,a,b) { r = (a+b); r= (r < p ? r : r-p); }
#define __GIVARO_ZPZIntType_N_ADD(r,p,a,b) { r = (a+b); if (r >= p) r-=p; }
// r += a
// #define __GIVARO_ZPZIntType_N_ADDIN(r,p,a) { r += a;  r= (r < p ? r : r-p); }
#define __GIVARO_ZPZIntType_N_ADDIN(r,p,a) { r += a;  if (r >= p) r-=p; }

// r <- a*b+c % p
// #define __GIVARO_ZPZIntType_N_MULADD(r,p,a,b,c) { r = (a*b+c) % p;  }
#define __GIVARO_ZPZIntType_N_MULADD(r,p,a,b,c) { r = a; r*=b; r+=c; r %= p;  }

// #define __GIVARO_ZPZIntType_N_MULADDIN(r,p,a,b) { r = (a*b+r) % p;  }
#define __GIVARO_ZPZIntType_N_MULADDIN(r,p,a,b) { r += (a*b); r %= p;  }

// a*b-c
//#define __GIVARO_ZPZIntType_N_MULSUB(r,p,a,b,c) { r = (a*b+p-c); r= (r<p ? r : r % p);  }
#define __GIVARO_ZPZIntType_N_MULSUB(r,p,a,b,c) { r = a*b; r+=p; r-=c; r= (r<p ? r : r % p);  }
// a*b-c
//#define __GIVARO_ZPZIntType_N_SUBMULIN(r,p,a,b) { r -= (a*b); if (r<0) { r+=p; r = (r<0 ? r % p : r); } }
#define __GIVARO_ZPZIntType_N_SUBMULIN(r,p,a,b) { r = p-r; r += a*b; r= (r<p ? r : r % p); __GIVARO_ZPZIntType_N_NEGIN(r,p); }

#define __GIVARO_ZPZIntType_N_NEG(r,p,a) { r = ( isZero(a) ? zero : p-a); }
#define __GIVARO_ZPZIntType_N_NEGIN(r,p) { r = ( isZero(r) ? zero : p-r); }

namespace Givaro
{
    // ------------------------- Arithmetic functions
    
    template<typename IntType, typename COMP> inline typename Modular<IntType, COMP>::Element&
    Modular<IntType, COMP>::mul (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_ZPZIntType_N_MUL(r,_p,a,b); return r;
    }

    template<typename IntType, typename COMP> inline typename Modular<IntType, COMP>::Element&
    Modular<IntType, COMP>::sub (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_ZPZIntType_N_SUB(r,_p,a,b); return r;
    }

    template<typename IntType, typename COMP> inline typename Modular<IntType, COMP>::Element&
    Modular<IntType, COMP>::add (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_ZPZIntType_N_ADD(r,_p,a,b); return r;
    }

    template<typename IntType, typename COMP> inline typename Modular<IntType, COMP>::Element&
    Modular<IntType, COMP>::neg (Element& r, const Element& a) const
    {
        __GIVARO_ZPZIntType_N_NEG(r,_p,a); return r;

    }

    template<typename IntType, typename COMP> inline typename Modular<IntType, COMP>::Element&
    Modular<IntType, COMP>::negin (Element& r) const
    {
        __GIVARO_ZPZIntType_N_NEGIN(r,_p);
        return r;
    }

    template<typename IntType, typename COMP> inline typename Modular<IntType, COMP>::Element&
    Modular<IntType, COMP>::inv (Element& u1, const Element& a) const
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

    template<typename IntType, typename COMP> inline typename Modular<IntType, COMP>::Element&
    Modular<IntType, COMP>::div (Element& r, const Element& a, const Element& b) const
    {
        typename Modular<IntType, COMP>::Element ib;
        inv(ib, b);
        __GIVARO_ZPZIntType_N_MUL(r,_p,a,ib);
        return r;
    }

    template<typename IntType, typename COMP> inline typename Modular<IntType, COMP>::Element&
    Modular<IntType, COMP>::mulin (Element& r, const Element& a) const
    {
        __GIVARO_ZPZIntType_N_MULIN(r,_p, a);
        return r;
    }

    template<typename IntType, typename COMP> inline typename Modular<IntType, COMP>::Element&
    Modular<IntType, COMP>::divin (Element& r, const Element& a) const
    {
        typename Modular<IntType, COMP>::Element ia;
        inv(ia, a);
        return mulin(r, ia);
    }

    template<typename IntType, typename COMP> inline typename Modular<IntType, COMP>::Element&
    Modular<IntType, COMP>::addin (Element& r, const Element& a) const
    {
        __GIVARO_ZPZIntType_N_ADDIN(r,_p, a);
        return r;
    }

    template<typename IntType, typename COMP> inline typename Modular<IntType, COMP>::Element&
    Modular<IntType, COMP>::subin (Element& r, const Element& a) const
    {
        __GIVARO_ZPZIntType_N_SUBIN(r,_p, a);
        return r;
    }


    template<typename IntType, typename COMP> inline typename Modular<IntType, COMP>::Element&
    Modular<IntType, COMP>::invin (Element& r) const
    {
        typename Modular<IntType, COMP>::Element t = r;
        return Modular<IntType, COMP>::inv(r,t);
    }

    template<typename IntType, typename COMP> inline typename Modular<IntType, COMP>::Element&
    Modular<IntType, COMP>::axpy (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        __GIVARO_ZPZIntType_N_MULADD(r, _p, a, b, c);
        return r;
    }

    template<typename IntType, typename COMP> inline typename Modular<IntType, COMP>::Element&
    Modular<IntType, COMP>::axpyin (Element& r, const Element& a, const Element& b) const
    {
        typename Modular<IntType, COMP>::Element tmp = r;
        __GIVARO_ZPZIntType_N_MULADDIN(tmp, _p, a, b);
        return r = (Element)tmp;
    }

    template<typename IntType, typename COMP> inline typename Modular<IntType, COMP>::Element&
    Modular<IntType, COMP>::axmy (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        __GIVARO_ZPZIntType_N_MULSUB(r, _p, a, b, c);
        return r;
    }

    // r = c - a*b
    template<typename IntType, typename COMP> inline typename Modular<IntType, COMP>::Element&
    Modular<IntType, COMP>::maxpy (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        typename Modular<IntType, COMP>::Element tmp = c;
        __GIVARO_ZPZIntType_N_SUBMULIN(tmp, _p, a, b );
        return r = (Element)tmp;
    }
    // r -= a*b
    template<typename IntType, typename COMP> inline typename Modular<IntType, COMP>::Element&
    Modular<IntType, COMP>::maxpyin (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_ZPZIntType_N_SUBMULIN(r, _p, a, b );
        return r;
    }
    // r = a*b - r
    template<typename IntType, typename COMP> inline typename Modular<IntType, COMP>::Element&
    Modular<IntType, COMP>::axmyin (Element& r, const Element& a, const Element& b) const
    {
        maxpyin(r,a,b);
        return negin(r);
    }

    //----- IO

    template<typename IntType, typename COMP> inline std::ostream&
    Modular<IntType, COMP>::write (std::ostream& s ) const
    {
        return s << "Modular<[IntType]> modulo " << residu();
    }

    template<typename IntType, typename COMP> inline std::istream&
    Modular<IntType, COMP>::read (std::istream& s, Element& a) const
    {
        s >> a;
        init(a, a);
        return s;
    }

    template<typename IntType, typename COMP> inline std::ostream&
    Modular<IntType, COMP>::write (std::ostream& s, const Element& a) const
    {
        return s << a;
    }

} // namespace Givaro

#endif
