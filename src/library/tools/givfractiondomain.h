// ===============================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <20 Jan 11 08:59:29 Jean-Guillaume.Dumas@imag.fr>
// Author: J-G. Dumas
// Description: fractions over a Ring domain
// ===============================================================
#ifndef __GIVARO_frac_domain_H
#define __GIVARO_frac_domain_H
#include <givaro/givpower.h>

#ifndef GIVABS
#define GIVABS(a) ((a)>0?(a):-(a))
#endif

namespace Givaro {
    template<class T1, class T2>
    struct Frac {
        T1 _num;
        T2 _den;
        Frac() {}
        template<class Z1, class Z2>
        Frac( const Z1& v1, const Z2& v2) : _num(v1), _den(v2) {}
        T1& nume() { return _num; }
        const T1& nume() const
        { return _num; }
        T2& deno() { return _den; }
        const T2& deno() const
        { return _den; }
    };


    template<class RingDom>
    struct FracDom : public RingDom {

    public :
        // -- Self_t
        typedef          FracDom<RingDom>		Self_t;

        // -- Exported types
        typedef	     RingDom			Ring_t;
        typedef typename RingDom::Element	Ring_E;
        typedef Frac<Ring_E, Ring_E>		Element;
        typedef Frac<Ring_E, Ring_E>		Rep;

        // -- Constantes
        const Rep zero;
        const Rep one;
        const Rep mOne;

        void reduce(Ring_E& a, Ring_E& b) const
        {
            Ring_E g; Ring_t::gcd(g,a,b);
            // Ring_t::write(std::cerr << "a BEF: ", a) << std::endl;
            // Ring_t::write(std::cerr << "b BEF: ", b) << std::endl;
            // Ring_t::write(std::cerr << "g GCD: ", g) << std::endl;
            if(! ( Ring_t::isOne(g) || Ring_t::isZero(g) ) ) {
                Ring_t::divin(a,g);
                Ring_t::divin(b,g);
            }
            // Ring_t::write(std::cerr << "a AFT: ", a) << std::endl;
            // Ring_t::write(std::cerr << "b AFT: ", b) << std::endl;
        }

        Rep& reduce(Rep& r) const
        {
            reduce(r._num,r._den);
            return r;
        }


        FracDom (const RingDom& R ) : Ring_t(R), zero(R.zero,R.one), one(R.one,R.one) , mOne(R.mOne,R.one){}
        FracDom (const Self_t& F) : Ring_t(static_cast<const Ring_t&>(F)), zero(F.zero), one(F.one), mOne(F.mOne) {}
        const Ring_t& getdomain() const
        { return static_cast<const Ring_t&>(*this); }
        const Ring_t& getring() const
        { return static_cast<const Ring_t&>(*this); }

        Rep& init(Rep& a) const
        { Ring_t::init(a._num); Ring_t::init(a._den); return a; }

        template<class XXX>
        Rep& init(Rep& p, const XXX &cste ) const
        {
            Ring_t::init(p._num,cste); Ring_t::assign(p._den,Ring_t::one);
            return p;
        }


        // -- Assignment p = q
        Rep& assign( Rep& p, const Rep& q) const
        {
            Ring_t::assign(p._num,q._num); Ring_t::assign(p._den,q._den);
            return p;
        }
        Rep& assign( Rep& p, const Ring_E& q) const
        {
            Ring_t::assign(p._num,q); Ring_t::assign(p._den,Ring_t::one);
            return p;
        }

        // -- Comparaison operator
        int isZero  ( const Rep& P ) const
        { return Ring_t::isZero(P._num); }
        int isOne   ( const Rep& P ) const
        { return Ring_t::areEqual(P._num, P._den); }
        int isMOne   ( const Rep& P ) const
        { return (Ring_t::areEqual(abs(P._num), abs(P._den)) && !isOne(P)); }

        int areEqual ( const Rep& P, const Rep& Q ) const
        {
            return Ring_t::areEqual(P._num, Q._num) && Ring_t::areEqual(P._den, Q._den) ; }
        int areNEqual( const Rep& P, const Rep& Q ) const
        {
            return Ring_t::areNEqual(P._num, Q._num) || Ring_t::areNEqual(P._den, Q._den) ;
        }

        // --
        std::istream& read ( std::istream& i ) {
            char tmp; i >> tmp; // '('
            Ring_t::read(i);
            return i >> tmp; // ')'
        }
        std::ostream& write( std::ostream& o ) const
        {
            return Ring_t::write(o << '(') << ')';
        }
        std::istream& read ( std::istream& i, Rep& n) const
        {
            char tmp;
            i >> tmp; GIVARO_ASSERT(tmp == '(', "Error in fraction read '(' not found");
            Ring_t::read(i,n._num);
            i >> tmp; GIVARO_ASSERT(tmp == ')', "Error in fraction read ')' not found");
            i >> tmp; GIVARO_ASSERT(tmp == '/', "Error in fraction read '/' not found");
            i >> tmp; GIVARO_ASSERT(tmp == '(', "Error in fraction read '(' not found");
            Ring_t::read(i,n._den);
            i >> tmp; GIVARO_ASSERT(tmp == ')', "Error in fraction read ')' not found");
            return i;
        }
        std::ostream& write( std::ostream& o, const Rep& n) const
        {
            return Ring_t::write(Ring_t::write(o << '(',n._num) << ")/(", n._den) << ')';
        }

        // -- Arithmetics operators
        Rep& mulin ( Rep& q, const Rep& a ) const
        {
            Ring_t::mulin(q._num,a._num);
            Ring_t::mulin(q._den,a._den);
            return reduce(q);
        }
        Rep& mulin ( Rep& q, const Ring_E& a ) const
        {
            Ring_E u(a);
            reduce(u,q._den);
            Ring_t::mulin(q._num,u);
            return q;
        }
        Rep& mul   ( Rep& q, const Rep& a, const Ring_E& b ) const
        {
            Ring_E u(b),v(a._den);
            reduce(u,v);
            Ring_t::mul(q._num,a._num,u);
            Ring_t::assign(q._den,v);
            return q;
        }

        Rep& mul   ( Rep& q, const Ring_E& a, const Rep& b ) const
        {
            return mul(q,b,a);
        }

        Rep& mul   ( Rep& q, const Rep& a, const Rep& b ) const
        {
            Ring_t::mul(q._num,a._num,b._num);
            Ring_t::mul(q._den,a._den,b._den);
            return reduce(q);
        }


        Rep& addin ( Rep& res, const Rep& u ) const
        {
            Ring_E g; Ring_t::gcd(g, res._den, u._den);
            Ring_E ud; Ring_t::div(ud, u._den, g);
            Ring_E vd; Ring_t::div(vd, res._den, g);
            Ring_t::mulin(res._num, ud);
            Ring_t::mulin(res._den, ud); // res *= uden/g
            Ring_t::mulin(vd, u._num);   // unum*= rden/g
            Ring_t::addin(res._num,vd);
            return reduce(res);
        }

        Rep& addin ( Rep& res, const Ring_E& a ) const
        {
            Ring_t::axpyin(res._num,res._den,a);
            return reduce(res);
        }

        Rep& add ( Rep& res, const Rep& u, const Rep& v ) const
        {
            Ring_E g; Ring_t::gcd(g, u._den, v._den);
            Ring_E ud; Ring_t::div(ud, u._den, g);
            Ring_E vd; Ring_t::div(vd, v._den, g);
            Ring_t::mul(res._num, u._num, vd);
            Ring_t::mul(res._den, u._den, vd); // res = u * vden/g
            Ring_t::mulin(ud, v._num);         // vnum* uden/g
            Ring_t::addin(res._num,ud);
            return reduce(res);
        }

        Rep& add ( Rep& res, const Rep& u, const Ring_E& a) const
        {
            Ring_t::axpy(res._num,u._num,u._den,a);
            Ring_t::assign(res._den,u._den);
            return reduce(res);
        }

        Rep& add ( Rep& res, const Ring_E& a, const Rep& u) const
        {
            return add(res,u,a);
        }

        Rep& subin ( Rep& res, const Rep& u ) const
        {
            Ring_E g; Ring_t::gcd(g, res._den, u._den);
            Ring_E ud; Ring_t::div(ud, u._den, g);
            Ring_E vd; Ring_t::div(vd, res._den, g);
            Ring_t::mulin(res._num, ud);
            Ring_t::mulin(res._den, ud); // res *= uden/g
            Ring_t::mulin(vd, u._num);   // unum*= rden/g
            Ring_t::subin(res._num,vd);
            return reduce(res);
        }

        Rep& subin ( Rep& res, const Ring_E& a ) const
        {
            Ring_t::maxpyin(res._num,res._den,a);
            return reduce(res);
        }

        Rep& sub ( Rep& res, const Rep& u, const Rep& v ) const
        {
            Ring_E g; Ring_t::gcd(g, u._den, v._den);
            Ring_E ud; Ring_t::div(ud, u._den, g);
            Ring_E vd; Ring_t::div(vd, v._den, g);
            Ring_t::mul(res._num, u._num, vd);
            Ring_t::mul(res._den, u._den, vd); // res = u * vden/g
            Ring_t::mulin(ud, v._num);         // vnum* uden/g
            Ring_t::subin(res._num,ud);
            return reduce(res);
        }

        Rep& sub ( Rep& res, const Rep& u, const Ring_E& a) const
        {
            Ring_t::maxpy(res._num,u._den,a,u._num);
            Ring_t::assign(res._den,u._den);
            return reduce(res);
        }

        Rep& sub ( Rep& res, const Ring_E& a, const Rep& u) const
        {
            Ring_t::axmy(res._num,a,u._den,u._num);
            Ring_t::assign(res._den,u._den);
            return reduce(res);
        }

        Rep& negin ( Rep& res ) const
        { Ring_t::negin(res._num); return res; }
        Rep& neg ( Rep& res, const Rep& u ) const
        {
            Ring_t::neg(res._num,u._num);
            Ring_t::assign(res._den,u._den);
            return res;
        }

        Rep& invin ( Rep& q) const
        {
            std::swap(q._num,q._den);
            return q;
        }

        Rep& inv( Rep& r, const Rep& u) const
        {
            Ring_t::assign(r._num, u._den);
            Ring_t::assign(r._den, u._num);
            return r;
        }
        Rep& inv(Rep& r, const Ring_E& a) {
            Ring_t::assign(r._den,a);
            Ring_t::assign(r._num,Ring_t::one);
            return r;
        }

        Rep& divin ( Rep& q, const Rep& a ) const
        {
            invin(q); mulin(q,a); return invin(q);
        }

        Rep& divin ( Rep& q, const Ring_E& a ) const
        {
            invin(q); mulin(q,a); return invin(q);
        }

        Rep& div   ( Rep& q, const Rep& a, const Rep& b ) const
        {
            inv(q,b); return mulin(q,a);
        }

        Rep& axpy  (Rep& r, const Rep& a, const Rep& x, const Rep& y) const
        {
            mul(r,x,y);
            return addin(r,a);
        }
        Rep& axpy  (Rep& r, const Rep& a, const Ring_E& x, const Rep& y) const
        {
            mul(r,x,y);
            return addin(r,a);
        }
        Rep& axpy  (Rep& r, const Rep& a, const Rep& x, const Ring_E& y) const
        {
            mul(r,x,y);
            return addin(r,a);
        }
        Rep& axpy  (Rep& r, const Ring_E& a, const Rep& x, const Rep& y) const
        {
            mul(r,x,y);
            return addin(r,a);
        }
        Rep& axpy  (Rep& r, const Ring_E& a, const Ring_E& x, const Rep& y) const
        {
            mul(r,x,y);
            return addin(r,a);
        }
        Rep& axpy  (Rep& r, const Ring_E& a, const Rep& x, const Ring_E& y) const
        {
            mul(r,x,y);
            return addin(r,a);
        }
        Rep& axpy  (Rep& r, const Rep& a, const Ring_E& x, const Ring_E& y) const
        {
            Ring_E m;Ring_t::mul(m,x,y);
            return add(r,a,m);
        }

        Rep& axpyin(Rep& r, const Rep& a, const Rep& x) const
        {
            Rep m; mul(m,a,x);
            return addin(r,m);
        }
        Rep& axpyin(Rep& r, const Ring_E& a, const Rep& x) const
        {
            Rep m; mul(m,a,x);
            return addin(r,m);
        }
        Rep& axpyin(Rep& r, const Rep& a, const Ring_E& x) const
        {
            Rep m; mul(m,a,x);
            return addin(r,m);
        }
        Rep& axpyin(Rep& r, const Ring_E& a, const Ring_E& x) const
        {
            Ring_E m; mul(m,a,x);
            return addin(r,m);
        }
        // -- maxpy: r <- y - a * x
        Rep& maxpy  (Rep& r, const Rep& a, const Rep& x, const Rep& y) const
        {
            mul(r,x,a); negin(r);
            return addin(r,y);
        }
        Rep& maxpy  (Rep& r, const Rep& a, const Ring_E& x, const Rep& y) const
        {
            mul(r,x,a); negin(r);
            return addin(r,y);
        }
        Rep& maxpy  (Rep& r, const Rep& a, const Rep& x, const Ring_E& y) const
        {
            mul(r,x,a); negin(r);
            return addin(r,y);
        }
        Rep& maxpy  (Rep& r, const Ring_E& a, const Rep& x, const Rep& y) const
        {
            mul(r,x,a); negin(r);
            return addin(r,y);
        }
        Rep& maxpy (Rep& r, const Ring_E& a, const Ring_E& x, const Rep& y) const
        {
            Ring_E m;Ring_t::mul(m,x,a);
            return sub(r,y,m);
        }
        Rep& maxpy  (Rep& r, const Ring_E& a, const Rep& x, const Ring_E& y) const
        {
            mul(r,x,a); negin(r);
            return addin(r,y);
        }
        Rep& maxpy  (Rep& r, const Rep& a, const Ring_E& x, const Ring_E& y) const
        {
            mul(r,x,a); negin(r);
            return addin(r,y);
        }

        // -- axmyin: r <-  a * x - r
        Rep& axmyin(Rep& r, const Rep& a, const Rep& x) const
        {
            maxpyin(r,a,x);
            return negin(r);
        }
        Rep& axmyin(Rep& r, const Ring_E& a, const Rep& x) const
        {
            maxpyin(r,a,x);
            return negin(r);
        }
        Rep& axmyin(Rep& r, const Rep& a, const Ring_E& x) const
        {
            maxpyin(r,a,x);
            return negin(r);
        }
        Rep& axmyin(Rep& r, const Ring_E& a, const Ring_E& x) const
        {
            maxpyin(r,a,x);
            return negin(r);
        }
        // -- maxpyin: r <- r - a * x
        Rep& maxpyin(Rep& r, const Rep& a, const Rep& x) const
        {
            Rep m; mul(m,a,x);
            return subin(r,m);
        }
        Rep& maxpyin(Rep& r, const Ring_E& a, const Rep& x) const
        {
            Rep m; mul(m,a,x);
            return subin(r,m);

        }
        Rep& maxpyin(Rep& r, const Rep& a, const Ring_E& x) const
        {
            Rep m; mul(m,a,x);
            return subin(r,m);

        }
        Rep& maxpyin(Rep& r, const Ring_E& a, const Ring_E& x) const
        {
            Ring_E m; mul(m,a,x);
            return subin(r,m);
        }
        // -- axmy: r <- a * x - y
        Rep& axmy  (Rep& r, const Rep& a, const Rep& x, const Rep& y) const
        {
            mul(r,a,x);
            return subin(r,y);
        }
        Rep& axmy  (Rep& r, const Ring_E& a, const Rep& x, const Rep& y) const
        {
            mul(r,a,x);
            return subin(r,y);
        }
        Rep& axmy  (Rep& r, const Rep& a, const Ring_E& x, const Rep& y) const
        {
            mul(r,a,x);
            return subin(r,y);
        }
        Rep& axmy  (Rep& r, const Rep& a, const Rep& x, const Ring_E& y) const
        {
            mul(r,a,x);
            return subin(r,y);
        }
        Rep& axmy  (Rep& r, const Ring_E& a, const Ring_E& x, const Rep& y) const
        {
            Ring_E m; mul(m,a,x);
            return sub(r,m,y);
        }
        Rep& axmy  (Rep& r, const Ring_E& a, const Rep& x, const Ring_E& y) const
        {
            mul(r,a,x);
            return subin(r,y);
        }
        Rep& axmy  (Rep& r, const Rep& a, const Ring_E& x, const Ring_E& y) const
        {
            mul(r,a,x);
            return subin(r,y);
        }

        // -- misc
        // -- W <-- P^n
        Rep& pow( Rep& W, const Rep& P, long n) const
        {
            unsigned long l = (unsigned long)GIVABS(n);
            if (n>0) {
                dom_power(W._num,P._num,l,static_cast<Ring_t&>(*this));
                dom_power(W._den,P._den,l,static_cast<Ring_t&>(*this));
            } else {
                dom_power(W._num,P._den,l,static_cast<Ring_t&>(*this));
                dom_power(W._den,P._num,l,static_cast<Ring_t&>(*this));
            }
            return W;
        }

        // -- Random generators
        template< class RandIter > Rep& random(RandIter& g, Rep& r) const
        {
            Ring_t::random(g, r._num);
            Ring_t::nonzerorandom(g, r._den);
            return r;
        }

        template< class RandIter, class XXX > Rep& random(RandIter& g, Rep& r, const XXX& s) const
        {
            Ring_t::random(g, r._num, s);
            Ring_t::nonzerorandom(g, r._den, s);
            return r;
        }

        template< class RandIter > Rep& nonzerorandom(RandIter& g, Rep& r) const
        {
            Ring_t::nonzerorandom(g, r._num);
            Ring_t::nonzerorandom(g, r._den);
            return r;
        }
        template< class RandIter, class XXX  > Rep& nonzerorandom(RandIter& g, Rep& r, const XXX& s) const
        {
            Ring_t::nonzerorandom(g, r._num, s);
            Ring_t::nonzerorandom(g, r._den, s);
            return r;
        }

    };

} // Givaro

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
