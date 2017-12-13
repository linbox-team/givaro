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
#ifndef __GIVARO_zpz_int_INL
#define __GIVARO_zpz_int_INL
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

    inline Modular<Integer>::Element&
    Modular<Integer>::mul (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_ZPZInteger_N_MUL(r,_p,a,b); return r;
    }

    inline Modular<Integer>::Element&
    Modular<Integer>::sub (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_ZPZInteger_N_SUB(r,_p,a,b); return r;
    }

    inline Modular<Integer>::Element&
    Modular<Integer>::add (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_ZPZInteger_N_ADD(r,_p,a,b); return r;
    }

    inline Modular<Integer>::Element&
    Modular<Integer>::neg (Element& r, const Element& a) const
    {
        __GIVARO_ZPZInteger_N_NEG(r,_p,a); return r;

    }

    inline Modular<Integer>::Element&
    Modular<Integer>::negin (Element& r) const
    {
        __GIVARO_ZPZInteger_N_NEGIN(r,_p);
        return r;
    }

    inline Modular<Integer>::Element&
    Modular<Integer>::inv (Element& r, const Element& a) const
    {
        return ::Givaro::inv(r,a,_p);
    }

    inline bool Modular<Integer>::isUnit(const Element& a) const 
    { 
        Element d;
        ::Givaro::gcd(d,a,_p); 
        return isOne(d) || isMOne(d); 
    }

    inline Modular<Integer>::Element&
    Modular<Integer>::div (Element& r, const Element& a, const Element& b) const
    {
        Element ib;
        inv(ib, b);
        __GIVARO_ZPZInteger_N_MUL(r,_p,a,ib);
        return r;
    }

    inline Modular<Integer>::Element&
    Modular<Integer>::mulin (Element& r, const Element& a) const
    {
        __GIVARO_ZPZInteger_N_MULIN(r,_p, a);
        return r;
    }

    inline Modular<Integer>::Element&
    Modular<Integer>::divin (Element& r, const Element& a) const
    {
        Modular<Integer>::Element ia;
        inv(ia, a);
        return mulin(r, ia);
    }

    inline Modular<Integer>::Element&
    Modular<Integer>::addin (Element& r, const Element& a) const
    {
        __GIVARO_ZPZInteger_N_ADDIN(r,_p, a);
        return r;
    }

    inline Modular<Integer>::Element&
    Modular<Integer>::subin (Element& r, const Element& a) const
    {
        __GIVARO_ZPZInteger_N_SUBIN(r,_p, a);
        return r;
    }

    inline Modular<Integer>::Element&
    Modular<Integer>::invin (Element& r) const
    {
        Element t = r;
        return ::Givaro::inv(r,t,_p);
    }

    inline Modular<Integer>::Element&
    Modular<Integer>::axpy (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        __GIVARO_ZPZInteger_N_MULADD(r, _p, a, b, c);
        return r;
    }

    inline Modular<Integer>::Element&
    Modular<Integer>::axpyin (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_ZPZInteger_N_MULADDIN(r, _p, a, b);
        return r;
    }

    inline Modular<Integer>::Element&
    Modular<Integer>::axmy (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        __GIVARO_ZPZInteger_N_MULSUB(r, _p, a, b, c);
        return r;
    }

    // r = c - a*b
    inline Modular<Integer>::Element&
    Modular<Integer>::maxpy (Element& r, const Element& a, const Element& b, const Element& c) const
    {
        Element tmp = c;
        __GIVARO_ZPZInteger_N_SUBMULIN(tmp, _p, a, b );
        return r = (Modular<Integer>::Element)tmp;
    }
    // r -= a*b
    inline Modular<Integer>::Element&
    Modular<Integer>::maxpyin (Element& r, const Element& a, const Element& b) const
    {
        __GIVARO_ZPZInteger_N_SUBMULIN(r, _p, a, b );
        return r;
    }
    // r = a*b - r
    inline Modular<Integer>::Element&
    Modular<Integer>::axmyin (Element& r, const Element& a, const Element& b) const
    {
        maxpyin(r,a,b);
        return negin(r);
    }

    // --------------------
    // ----- Initialisation
    
    inline typename Modular<Integer>::Element&
    Modular<Integer>::init(Element& x) const
    {
        return x = zero;
    }

    inline typename Modular<Integer>::Element&
    Modular<Integer>::assign ( Element& r, const Element& a ) const
    {
	return r = a;
    }

    // ------------
    // ----- Reduce

    inline typename Modular<Integer>::Element&
    Modular<Integer>::reduce(Element& r, const Element& a) const
    {
	r = a % static_cast<Element>(_p);
	if (r < 0) r = static_cast<Element>(r + _p);
	return r;
    }

    inline typename Modular<Integer>::Element&
    Modular<Integer>::reduce(Element& r) const
    {
	r %= static_cast<Element>(_p);
	if (r < 0) r = static_cast<Element>(r + _p);
	return r;
    }

    //----- IO
    
    inline std::ostream& Modular<Integer>::write (std::ostream& s ) const
    {
        return s << "Modular<Integer> modulo " << residu();
    }

    inline std::istream& Modular<Integer>::read (std::istream& s, Element& a) const
    {
        s >> a;
        init(a, a);
        return s;
    }

    inline std::ostream& Modular<Integer>::write (std::ostream& s, const Element& a) const
    {
        return s << a;
    }

} // namespace Givaro

#endif // __GIVARO_zpz_int_INL

