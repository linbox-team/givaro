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

#define __GIVARO_ZPZIntType_N_NEG(r,p,a) { r = ( Modular<IntType, COMP, Enable>::isZero(a) ? zero : p-a); }
#define __GIVARO_ZPZIntType_N_NEGIN(r,p) { r = ( isZero(r) ? zero : p-r); }

namespace Givaro
{
    // ------------------------- Arithmetic functions
    
    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::mul (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_ZPZIntType_N_MUL(r,_p,a,b); return r;
    }

    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::sub (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_ZPZIntType_N_SUB(r,_p,a,b); return r;
    }

    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::add (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_ZPZIntType_N_ADD(r,_p,a,b); return r;
    }

    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::neg (Element& r, const Element& a) const
    {
        //__GIVARO_ZPZIntType_N_NEG(r,_p,a); return r;
        return r = (a == zero ? zero : _p - a);
    }

    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::negin (Element& r) const
    {
        __GIVARO_ZPZIntType_N_NEGIN(r,_p);
        return r;
    }

    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::inv (Element& u1, const Element& a) const
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

    
//    template<typename IntType, typename COMP, typename Enable> 
//    inline bool Modular<IntType, COMP, Enable>::isUnit(const Element& a) const 
//    { 
//        Element u,d; 
//        invext(u,d,a,_p); 
//        return isOne(d) || isMOne(d); 
//    }


    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::div (Element& r, const Element& a, const Element& b) const
    {
        typename Modular<IntType, COMP, Enable>::Element ib;
        inv(ib, b);
        __GIVARO_ZPZIntType_N_MUL(r,_p,a,ib);
        return r;
    }

    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::mulin (Element& r, const Element& a) const
    {
        __GIVARO_ZPZIntType_N_MULIN(r,_p, a);
        return r;
    }

    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::divin (Element& r, const Element& a) const
    {
        typename Modular<IntType, COMP, Enable>::Element ia;
        inv(ia, a);
        return mulin(r, ia);
    }

    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::addin (Element& r, const Element& a) const
    {
        __GIVARO_ZPZIntType_N_ADDIN(r,_p, a);
        return r;
    }

    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::subin (Element& r, const Element& a) const
    {
        __GIVARO_ZPZIntType_N_SUBIN(r,_p, a);
        return r;
    }


    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::invin (Element& r) const
    {
        typename Modular<IntType, COMP, Enable>::Element t = r;
        return Modular<IntType, COMP, Enable>::inv(r,t);
    }

    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::axpy (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        __GIVARO_ZPZIntType_N_MULADD(r, _p, a, b, c);
        return r;
    }

    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::axpyin (Element& r, const Element& a, const Element& b) const
    {
        typename Modular<IntType, COMP, Enable>::Element tmp = r;
        __GIVARO_ZPZIntType_N_MULADDIN(tmp, _p, a, b);
        return r = (Element)tmp;
    }

    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::axmy (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        __GIVARO_ZPZIntType_N_MULSUB(r, _p, a, b, c);
        return r;
    }

    // r = c - a*b
    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::maxpy (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        typename Modular<IntType, COMP, Enable>::Element tmp = c;
        __GIVARO_ZPZIntType_N_SUBMULIN(tmp, _p, a, b );
        return r = (Element)tmp;
    }
    // r -= a*b
    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::maxpyin (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_ZPZIntType_N_SUBMULIN(r, _p, a, b );
        return r;
    }
    // r = a*b - r
    template<typename IntType, typename COMP, typename Enable> inline typename Modular<IntType, COMP, Enable>::Element&
    Modular<IntType, COMP, Enable>::axmyin (Element& r, const Element& a, const Element& b) const
    {
        maxpyin(r,a,b);
        return negin(r);
    }


} // namespace Givaro

#endif
