// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givzpzInt.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: JG Dumas
// $Id: givzpzInt.inl,v 1.11 2011-02-02 17:16:43 briceboyer Exp $
// ==========================================================================
#ifndef __GIVARO_mod_integer_INL
#define __GIVARO_mod_integer_INL
// Description:

// ---------
// -- normalized operations
// ---------

// r = a*b
#define __GIVARO_ZPZInteger_N_MUL(r,p,a,b) { Integer::mul(r,a,b); Integer::modin(r,p); }
// r *= a
#define __GIVARO_ZPZInteger_N_MULIN(r,p,a) {  Integer::mulin(r,a); Integer::modin(r,p);  }

// r = a - b
#define __GIVARO_ZPZInteger_N_SUB(r,p,a,b) { Integer::sub(r,a,b); if (sign(r) < 0) Integer::addin(r,p); }
// r -= a
#define __GIVARO_ZPZInteger_N_SUBIN(r,p,a) { Integer::subin(r,a) ; if ( sign(r) < 0) Integer::addin(r,p); }

// r = a+b
#define __GIVARO_ZPZInteger_N_ADD(r,p,a,b) { Integer::add(r,a,b); if (r >= p) Integer::subin(r,p); }
// r += a
#define __GIVARO_ZPZInteger_N_ADDIN(r,p,a) { Integer::addin(r,a);  if (r >= p) Integer::subin(r,p); }

// r <- a*b+c % p
#define __GIVARO_ZPZInteger_N_MULADD(r,p,a,b,c) { Integer::axpy(r,a,b,c); Integer::modin(r,p);  }
#define __GIVARO_ZPZInteger_N_MULADDIN(r,p,a,b) { Integer::axpyin(r,a,b); Integer::modin(r,p); }

// a*b-c
#define __GIVARO_ZPZInteger_N_MULSUB(r,p,a,b,c) { Integer::axmy(r,a,b,c); Integer::modin(r,p);  }
// a*b-c
#define __GIVARO_ZPZInteger_N_SUBMULIN(r,p,a,b) { Integer::maxpyin(r,a,b); Integer::modin(r,p) ; }

#define __GIVARO_ZPZInteger_N_NEG(r,p,a) { if (isZero(a)) r=a; else Integer::sub(r,p,a);  }
#define __GIVARO_ZPZInteger_N_NEGIN(r,p) { if (! isZero(r)) Integer::sub(r,p,r); }

namespace Givaro {
    
    // ------------------------- Arithmetic functions

    inline typename Mod<Integer>::Element&
    Mod<Integer>::mul (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_ZPZInteger_N_MUL(r,_p,a,b); return r;
    }

    inline typename Mod<Integer>::Element&
    Mod<Integer>::sub (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_ZPZInteger_N_SUB(r,_p,a,b); return r;
    }

    inline typename Mod<Integer>::Element&
    Mod<Integer>::add (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_ZPZInteger_N_ADD(r,_p,a,b); return r;
    }

    inline typename Mod<Integer>::Element&
    Mod<Integer>::neg (Element& r, const Element& a) const
    {
        __GIVARO_ZPZInteger_N_NEG(r,_p,a); return r;

    }

    inline typename Mod<Integer>::Element&
    Mod<Integer>::negin (Element& r) const
    {
        __GIVARO_ZPZInteger_N_NEGIN(r,_p);
        return r;
    }

    inline typename Mod<Integer>::Element&
    Mod<Integer>::inv (Element& r, const Element& a) const
    {
        return ::Givaro::inv(r,a,_p);
    }

//    TMPL
//    inline bool this->Mod<Integer>::isUnit(const Element& a) const 
//    { 
//        Element d;
//        ::Givaro::gcd(d,a,_p); 
//        return isOne(d) || isMOne(d); 
//    }

    inline typename Mod<Integer>::Element&
    Mod<Integer>::div (Element& r, const Element& a, const Element& b) const
    {
        Element ib;
        inv(ib, b);
        __GIVARO_ZPZInteger_N_MUL(r,_p,a,ib);
        return r;
    }

    inline typename Mod<Integer>::Element&
    Mod<Integer>::mulin (Element& r, const Element& a) const
    {
        __GIVARO_ZPZInteger_N_MULIN(r,_p, a);
        return r;
    }

    inline typename Mod<Integer>::Element&
    Mod<Integer>::divin (Element& r, const Element& a) const
    {
        Mod<Integer>::Element ia;
        inv(ia, a);
        return mulin(r, ia);
    }

    inline typename Mod<Integer>::Element&
    Mod<Integer>::addin (Element& r, const Element& a) const
    {
        __GIVARO_ZPZInteger_N_ADDIN(r,_p, a);
        return r;
    }

    inline typename Mod<Integer>::Element&
    Mod<Integer>::subin (Element& r, const Element& a) const
    {
        __GIVARO_ZPZInteger_N_SUBIN(r,_p, a);
        return r;
    }

    inline typename Mod<Integer>::Element&
    Mod<Integer>::invin (Element& r) const
    {
        Element t = r;
        return ::Givaro::inv(r,t,_p);
    }

    inline typename Mod<Integer>::Element&
    Mod<Integer>::axpy (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        __GIVARO_ZPZInteger_N_MULADD(r, _p, a, b, c);
        return r;
    }

    inline typename Mod<Integer>::Element&
    Mod<Integer>::axpyin (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_ZPZInteger_N_MULADDIN(r, _p, a, b);
        return r;
    }

    inline typename Mod<Integer>::Element&
    Mod<Integer>::axmy (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        __GIVARO_ZPZInteger_N_MULSUB(r, _p, a, b, c);
        return r;
    }

    // r = c - a*b
    inline typename Mod<Integer>::Element&
    Mod<Integer>::maxpy (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        Element tmp = c;
        __GIVARO_ZPZInteger_N_SUBMULIN(tmp, _p, a, b );
        return r = (Mod<Integer>::Element)tmp;
    }
    // r -= a*b
    inline typename Mod<Integer>::Element&
    Mod<Integer>::maxpyin (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_ZPZInteger_N_SUBMULIN(r, _p, a, b );
        return r;
    }
    // r = a*b - r
    inline typename Mod<Integer>::Element&
    Mod<Integer>::axmyin (Element& r, const Element& a, const Element& b) const
    {
        maxpyin(r,a,b);
        return negin(r);
    }

    // --------------------
    // ----- Initialisation
    
    inline typename Mod<Integer>::Element&
    Mod<Integer>::init(Element& x) const
    {
        return x = this->zero;
    }

    // ------------
    // ----- Reduce

    inline typename Mod<Integer>::Element&
    Mod<Integer>::reduce(Element& r, const Element& a) const
    {
	r = a % _p;
	if (r < 0) r = r + _p;
	return r;
    }

    inline typename Mod<Integer>::Element&
    Mod<Integer>::reduce(Element& r) const
    {
	r %= _p;
	if (r < 0) r = r + _p;
	return r;
    }


} // namespace Givaro

#endif // __GIVARO_mod_integer_INL

