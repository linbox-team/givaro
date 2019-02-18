// ===============================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <19 Oct 10 18:35:20 Jean-Guillaume.Dumas@imag.fr>
// Author: J-G. Dumas
// Description: Quotients over a Ring domain
// ===============================================================
#ifndef __GIVARO_quotient_domain_H
#define __GIVARO_quotient_domain_H
#include <givaro/givpower.h>
#ifndef GIVABS
#define GIVABS(a) ((a)>0?(a):-(a))
#endif

namespace Givaro {
    template<class RingDom>
    struct QuotientDom : public RingDom {
    public :
        // -- Self_t
        typedef          QuotientDom<RingDom>	Self_t;

        // -- Exported types
        typedef	     RingDom			Ring_t;
        typedef typename RingDom::Element		Ring_E;
        typedef Ring_E				Element;
        typedef Ring_E				Rep;
    protected :
        Rep _modulo;

    public :
        QuotientDom (const RingDom& R, const Element& Mod ) : Ring_t(R), _modulo(Mod) {}
        QuotientDom (const Self_t& F) : Ring_t(static_cast<const Ring_t&>(F)), _modulo(F._modulo) {}

        Rep& init(Rep& a) const
        { return Ring_t::modin(Ring_t::init(a),_modulo); }

        template<class XXX>
        Rep& init(Rep& p, const XXX &cste ) const
        {
            return Ring_t::modin(Ring_t::init(p,cste),_modulo);
        }

        Rep& assign(Rep& p) const
        {
            return Ring_t::modin(p,_modulo);
        }
        Rep& assign(Rep& p, const Rep& Q) const
        {
            return Ring_t::modin(Ring_t::assign(p,Q),_modulo);
        }

        // -- Comparaison operator
        int isZero  ( const Rep& P ) const
        { return Ring_t::isZero(P); }
        int isOne   ( const Rep& P ) const
        { return Ring_t::isOne(P); }
        int isMOne   ( const Rep& P ) const
        { return Ring_t::isMOne(P); }

        int areEqual ( const Rep& P, const Rep& Q ) const
        {
            return Ring_t::areEqual(P, Q);
        }
        int areNEqual( const Rep& P, const Rep& Q ) const
        {
            return Ring_t::areNEqual(P, Q) ;
        }

        // --
        std::istream& read ( std::istream& i ) {
            char tmp;
            return Ring_t::read(Ring_t::read(i) >> tmp);
        }
        std::ostream& write( std::ostream& o ) const
        {
            return Ring_t::write( Ring_t::write(o) << '/', _modulo);
        }
        std::istream& read ( std::istream& i, Rep& n) const
        {
            return Ring_t::read(i,n);
        }
        std::ostream& write( std::ostream& o, const Rep& n) const
        {
            return Ring_t::write(o,n);
        }

        // -- Arithmetics operators
        Rep& mulin ( Rep& q, const Rep& a ) const
        {
            return Ring_t::modin(Ring_t::mulin(q,a), _modulo);
        }
        Rep& mul   ( Rep& q, const Rep& a, const Rep& b ) const
        {
            return Ring_t::modin(Ring_t::mul(q,a,b), _modulo);
        }
        Rep& addin ( Rep& r, const Rep& u ) const
        {
            return Ring_t::modin(Ring_t::addin(r,u), _modulo);
        }
        Rep& add ( Rep& r, const Rep& u, const Rep& v ) const
        {
            return Ring_t::modin(Ring_t::add(r,u,v), _modulo);
        }
        Rep& subin ( Rep& r, const Rep& u ) const
        {
            return Ring_t::modin(Ring_t::subin(r,u), _modulo);
        }
        Rep& sub ( Rep& r, const Rep& u, const Rep& v ) const
        {
            return Ring_t::modin(Ring_t::sub(r,u,v), _modulo);
        }
        Rep& negin ( Rep& r ) const
        {
            return Ring_t::modin(Ring_t::negin(r),_modulo);
        }
        Rep& neg ( Rep& r, const Rep& u ) const
        {
            return Ring_t::modin(Ring_t::neg(r,u),_modulo);
        }
        Rep& invin ( Rep& q) const
        {
            Rep t; Ring_t::invmod(t,q,_modulo);
            return Ring_t::assign(q,t);
        }
        Rep& inv( Rep& r, const Rep& u) const
        {
            return Ring_t::invmod(r,u,_modulo);
        }

        Rep& divin ( Rep& q, const Rep& a ) const
        {
            Rep t;
            return this->mulin(q,this->inv(t,a));
        }
        Rep& div   ( Rep& q, const Rep& a, const Rep& b ) const
        {
            return this->mulin(this->inv(q, b),a);
        }
        Rep& axpy  (Rep& r, const Rep& a, const Rep& x, const Rep& y) const
        {
            return Ring_t::modin(Ring_t::axpy(r,a,x,y), _modulo);
        }
        Rep& axpyin(Rep& r, const Rep& a, const Rep& x) const
        {
            return Ring_t::modin(Ring_t::axpyin(r,a,x), _modulo);
        }
        // -- maxpy: r <- y - a * x
        Rep& maxpy  (Rep& r, const Rep& a, const Rep& x, const Rep& y) const
        {
            return Ring_t::modin(Ring_t::maxpy(r,a,x,y), _modulo);
        }
        // -- axmyin: r <- a * x - r
        Rep& axmyin(Rep& r, const Rep& a, const Rep& x) const
        {
            return Ring_t::modin(Ring_t::axmyin(r,a,x), _modulo);
        }
        // -- maxpyin: r <- r - a * x
        Rep& maxpyin(Rep& r, const Rep& a, const Rep& x) const
        {
            return Ring_t::modin(Ring_t::maxpyin(r,a,x), _modulo);
        }
        // -- axmy: r <- a * x - y
        Rep& axmy  (Rep& r, const Rep& a, const Rep& x, const Rep& y) const
        {
            return Ring_t::modin(Ring_t::axmy(r,a,x,y), _modulo);
        }
        // -- misc
        // -- W <-- P^n
        Rep& pow( Rep& W, const Rep& P, long n) const
        {
            unsigned long l = (unsigned long)GIVABS(n);
            if (n>0)
                return dom_power(W, P, l, *this);
            else {
                Rep invP; this->inv(invP,P);
                return dom_power(W, invP, l, *this);
            }
        }
        // -- Random generators
        template< class RandIter >
        Rep& random(RandIter& g, Rep& r) const
        {
            return Ring_t::modin(Ring_t::random(g, r),_modulo);
        }

        template< class RandIter, class XXX >
        Rep& random(RandIter& g, Rep& r, const XXX& s) const
        {
            return Ring_t::modin(Ring_t::random(g, r, s),_modulo);
        }

        template< class RandIter > Rep&
        nonzerorandom(RandIter& g, Rep& r) const
        {
            return Ring_t::modin(Ring_t::nonzerorandom(g, r),_modulo);
        }
        template< class RandIter, class XXX  >
        Rep& nonzerorandom(RandIter& g, Rep& r, const XXX& s) const
        {
            return Ring_t::modin(Ring_t::nonzerorandom(g, r, s),_modulo);
        }

    };

} // Givaro

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
