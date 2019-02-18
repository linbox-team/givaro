// ===============================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <18 Feb 11 16:04:19 Jean-Guillaume.Dumas@imag.fr>
// Author: J-G. Dumas
// Description: Pieces of polynomials as defined in
// [Arne Storjohann, High-Order Lifting
//  ISSAC'2002, pp. 246-254, ACM Press, July 2002].
// ===============================================================
#ifndef __GIVARO_trunc_domain_H
#define __GIVARO_trunc_domain_H
#include <givaro/givpoly1dense.h>
#ifndef __PATHCC__
#include <utility>
#endif


namespace Givaro {

    template <class Domain>
    class TruncDom : public Poly1Dom<Domain,Dense> {
    public :
        // -- Self_t
        typedef          TruncDom<Domain>	Self_t;
        // -- Father_t
        typedef          Poly1Dom<Domain,Dense>  	Father_t;
        // -- Exported types
        typedef          Domain		       	Domain_t;
        typedef typename Domain::Element      Type_t;
        typedef typename Father_t::Storage_t	Polynomial_t;
        typedef std::pair<Polynomial_t, Degree>	Storage_t;
        typedef          Storage_t                Rep;
        typedef          Storage_t                Element;


        Storage_t zero, one,mOne;

        TruncDom (const Domain& d, const Indeter& X = Indeter() ) : Father_t(d,X) {
            this->assign(zero,Father_t::zero);
            this->assign(one,Father_t::one);
            this->assign(mOne,Father_t::mOne);
        }
        TruncDom (const Self_t& t) : Father_t(static_cast<const Father_t&>(t)) {
            this->assign(zero,Father_t::zero);
            this->assign(one,Father_t::one);
            this->assign(mOne,Father_t::mOne);
        }
        TruncDom (const Father_t& t) : Father_t(t) {
            this->assign(zero,Father_t::zero);
            this->assign(one,Father_t::one);
            this->assign(mOne,Father_t::mOne);
        }

        const Father_t& getpoldomain() const
        {
            return static_cast<const Father_t&>(*this);
        }


        Rep& init(Rep& p) const
        { Father_t::init(p.first); p.second=0; return p; }

        template<class XXX>
        Rep& init(Rep& p, const XXX &cste ) const
        {
            Father_t::init(p.first,cste); p.second=0; return p;
        }

        // -- For polynomial = lcoeff X^deg
        template<class XXX>
        Rep& init (Rep& p, const Degree deg , const XXX& lcoeff) const
        {
            Father_t::init(p.first,Degree(0),lcoeff); p.second=deg; return p;
        }

        Polynomial_t& convert (Polynomial_t& r, const Rep& P) const
        {
            Father_t::assign(r, P.first);
            Polynomial_t mon;
            Father_t::init(mon, P.second);
            return Father_t::mulin(r, mon);
        }

        template<class XXX>
        XXX& convert (XXX& r, const Rep& P) const
        {
            Polynomial_t eP; Father_t::init(eP);
            this->convert(eP, P);
            return Father_t::convert(r, eP);
        }

        //    F.assign(P[deg], lcoeff);
        Rep& assign (Rep& p, const Degree deg , const Type_t& lcoeff) const
        {
            Father_t::assign(p.first,Degree(0),lcoeff);
            Father_t::setdegree(p.first);
            p.second=deg;
            return setval(p);
        }

        // -- Assign polynomial with field value : F.assign(p[0],cste)
        Rep& assign(Rep& p, const Type_t &cste ) const
        {
            return assign(p, Degree(0), cste);
        }
        // -- Assignment p = q
        Rep& assign( Rep& p, const Rep& q) const
        {
            Father_t::assign(p.first, q.first);
            p.second = q.second;
            return p;
        }

        Rep& assign(Rep& p, const Polynomial_t& r ) const
        {
            Father_t::assign(p.first,r); p.second=0; return setval(p);
        }

        Rep& assign(Rep& p, const Polynomial_t& r, const Degree v, const Degree d) const
        {
            Father_t::assign(p.first,r); p.second=0; return truncin(p,v,d);
        }


        Rep& mulin(Rep& p, const Degree& s) const
        {
            p.second += s; return p;
        }

        Rep& divin(Rep& p, const Degree& s) const
        {
            p.second -= s; p.second=(p.second<0?0:p.second); return p;
        }

        // -- Compute the degree of P
        Rep& setdegree( Rep& P ) const
        {
            Father_t::setdegree(P.first);
            if (P.first.size() <= 0) P.second=0;
            return P;
        }

        // -- Compute the valuation of P
        Rep& setval( Rep& P ) const
        {
            setdegree(P);
            if (P.first.size() <= 0) return P;
            typename Polynomial_t::iterator it = P.first.begin();
            if (! this->_domain.isZero(*it)) return P;
            for(++it,++P.second; it != P.first.end(); ++it,++P.second) {
                if (! this->_domain.isZero(*it)) {
                    P.first.erase(P.first.begin(),it);
                    return P;
                }
            }
            P.first.resize(0); P.second=0;
            return P;
        }

        // -- Returns the degree of polynomial
        Degree& degree(Degree& d, const Rep& P) const
        {
            Father_t::degree(d, P.first);
            return d+=P.second;
        }

        // -- Returns the valuation of polynomial
        Degree& val(Degree& d, const Rep& P) const
        {
            this->setval(const_cast<Rep&>(P));
            return d=P.second;
        }


        // -- Comparaison operator
        int isZero  ( const Rep& P ) const
        { return Father_t::isZero(P.first); }
        int isOne   ( const Rep& P ) const
        {
            Degree vP;val(vP,P);
            return ((vP == 0) && Father_t::isOne(P.first));
        }
        int isMOne   ( const Rep& P ) const
        {
            Degree vP;val(vP,P);
            return ((vP == 0) && Father_t::isMOne(P.first));
        }


        int areEqual ( const Rep& P, const Rep& Q ) const
        {
            Degree vP;val(vP,P);
            Degree vQ;val(vQ,Q);
            return (vP == vQ) &&  Father_t::areEqual(P.first,Q.first);
        }

        int areNEqual( const Rep& P, const Rep& Q ) const
        {
            Degree vP;val(vP,P);
            Degree vQ;val(vQ,Q);
            return (vP != vQ) ||  Father_t::areNEqual(P.first,Q.first);
        }


        Rep& shift(Rep& p, const Degree& s) const;
        Rep& truncin(Rep& p, const Degree& v, const Degree& d) const;
        Rep& trunc(Rep& p, const Rep& R, const Degree& v, const Degree& d) const
        {
            return this->truncin(this->assign(p,R),v,d);
        }

        // I/O
        std::istream& read ( std::istream& i ) {
            char tmp, t[5];
            return Father_t::read(i>> tmp)>> t ;
        }
        std::ostream& write( std::ostream& o ) const
        {
            return Father_t::write(o << '[') << "]_i^j";
        }
        std::istream& read ( std::istream& i, Rep& n) const
        {
            char tmp;
            return Father_t::read(i>>tmp,n.first)>>tmp>>tmp>>tmp>>tmp>> n.second;
        }
        std::ostream& write( std::ostream& o, const Rep& n) const
        {

            return Father_t::write(o<<'(',n.first)<<")*" << this->_x << '^' << n.second;
        }

        Rep& expand(Rep& P, const Degree& d) const
        {
            Degree vP; val(vP, P);
            if (vP > d) {
                P.first.insert(P.first.begin(),(size_t)value(vP-d),this->_domain.zero);
                P.second = d;
            }
            return P;
        }


        // -- Arithmetics operators
        Rep& addin ( Rep& R, const Rep& P) const;
        Rep& add ( Rep& res, const Rep& u, const Rep& v ) const
        {
            assign(res,u);
            return addin(res,v);
        }

        Rep& addin ( Rep& R, const Rep& P, const Degree& v, const Degree& d) const;
        Rep& add ( Rep& res, const Rep& u, const Rep& v, const Degree& Val, const Degree& deg) const
        {
            assign(res,u);
            return addin(res,v,Val,deg);
        }

        Rep& neg(Rep& R, const Rep& P) const
        {
            Father_t::neg(R.first,P.first);
            R.second = P.second;
            return R;
        }
        Rep& negin(Rep& R) const
        {
            Father_t::negin(R.first);
            return R;
        }

        Rep& sub ( Rep& res, const Rep& u, const Rep& v ) const
        {
            assign(res,u);
            return this->subin(res,v);
        }
        Rep& subin ( Rep& R, const Rep& P) const ;
        Rep& sub ( Rep& R, const Rep& P, const Rep& Q, const Degree& v, const Degree& d) const
        {
            return this->addin(this->neg(R,Q),P,v,d);

        }

        Rep& subin ( Rep& R, const Rep& P, const Degree& v, const Degree& d) const
        {
            return this->negin(this->addin(this->negin(R),P,v,d));
        }
        Rep& mul ( Rep& res, const Rep& u, const Rep& v ) const
        {
            Father_t::mul(res.first,u.first,v.first);
            res.second = u.second+v.second;
            return res;
        }
        Rep& mulin ( Rep& P, const Rep& Q ) const
        {
            Father_t::mulin(P.first,Q.first);
            P.second += Q.second;
            return P;
        }


        Rep& mul( Rep& r, const Rep& u, const Rep& v, const Degree& Val, const Degree& deg) const;
        Rep& mulin( Rep& r, const Rep& v, const Degree& Val, const Degree& deg) const
        {
            Rep tmp(r);
            return mul(r,tmp,v,Val,deg);
        }

        Rep& axpy  (Rep& r, const Rep& a, const Rep& x, const Rep& y) const
        {
            return this->addin( this->mul(r,a,x), y);
        }
        Rep& axpyin(Rep& r, const Rep& a, const Rep& x) const
        {
            Rep tmp; this->init(tmp);
            return this->addin(r, this->mul(tmp,a,x));
        }
        Rep& axpy  (Rep& r, const Rep& a, const Rep& x, const Rep& y, const Degree& Val, const Degree& deg) const
        {
            return this->addin(this->mul(r,a,x,Val,deg), y, Val, deg);
        }
        Rep& axpyin  (Rep& r, const Rep& a, const Rep& x, const Degree& Val, const Degree& deg) const
        {
            Rep tmp; this->init(tmp);
            return this->addin(r, this->mul(tmp,a,x,Val,deg), Val, deg);
        }

        Rep& axmy  (Rep& r, const Rep& a, const Rep& x, const Rep& y) const
        {
            return this->subin( this->mul(r,a,x), y);
        }
        Rep& axmyin(Rep& r, const Rep& a, const Rep& x) const
        {
            return this->negin(this->maxpyin(r,a,x));
        }
        Rep& axmy  (Rep& r, const Rep& a, const Rep& x, const Rep& y, const Degree& Val, const Degree& deg) const
        {
            return this->subin(this->mul(r,a,x,Val,deg), y, Val, deg);
        }
        Rep& axmyin  (Rep& r, const Rep& a, const Rep& x, const Degree& Val, const Degree& deg) const
        {
            return this->negin(this->maxpyin(r,a,x,Val,deg));
        }

        Rep& maxpy  (Rep& r, const Rep& a, const Rep& x, const Rep& y) const
        {
            return this->addin( this->negin(this->mul(r,a,x)), y);
        }
        Rep& maxpyin(Rep& r, const Rep& a, const Rep& x) const
        {
            Rep tmp;
            return this->subin(r, this->mul(tmp,a,x));
        }
        Rep& maxpy  (Rep& r, const Rep& a, const Rep& x, const Rep& y, const Degree& Val, const Degree& deg) const
        {
            return this->addin( this->negin(this->mul(r,a,x,Val,deg)), y, Val, deg);
        }
        Rep& maxpyin  (Rep& r, const Rep& a, const Rep& x, const Degree& Val, const Degree& deg) const
        {
            Rep tmp;
            return this->subin(r, this->mul(tmp,a,x,Val,deg), Val, deg);
        }

        // -- Random generators
        // -- Random dense polynomial of degree 0
        template< class RandIter > Rep& random(RandIter& g, Rep& r) const
        {
            Father_t::random(g,r.first);
            r.second = rand();
            return r;
        }

        // -- Random dense polynomial of size s
        template< class RandIter > Rep& random(RandIter& g, Rep& r, long s) const
        {
            Father_t::random(g,r.first,s);
            r.second = rand() % s;
            return r;
        }
        // -- Random dense polynomial of degree d
        template< class RandIter > Rep& random(RandIter& g, Rep& r, Degree s) const
        {
            Father_t::random(g,r.first,s);
            r.second = rand() % s.value();
            return r;
        }

        Rep& random(GivRandom& g, Rep& r, Degree s) const
        {
            Father_t::random(g,r.first,s);
            r.second = (Degree)(long)((unsigned long)g() % (unsigned long)((s.value()<<1)|1));
            return r;
        }
        // -- Random dense polynomial with same size as b.
        template< class RandIter > Rep& random(RandIter& g, Rep& r, const Rep& b) const
        {
            Father_t::random(g,r.first,b);
            r.second = b.second;
            return r;
        }

        template< class RandIter > Rep& nonzerorandom(RandIter& g, Rep& r) const{
            Father_t::nonzerorandom(g,r.first);
            r.second = rand();
            return r;
        }
        template< class RandIter > Rep& nonzerorandom(RandIter& g, Rep& r, long s) const
        {
            Father_t::nonzerorandom(g,r.first,s);
            r.second = rand() % s;
            return r;
        }
        template< class RandIter > Rep& nonzerorandom(RandIter& g, Rep& r, Degree s) const
        {
            Father_t::nonzerorandom(g,r.first,s);
            r.second = rand() % s.value();
            return r;
        }
        template< class RandIter > Rep& nonzerorandom(RandIter& g, Rep& r, const Rep& b) const
        {
            Father_t::random(g,r.first,b);
            r.second = b.second;
            return r;
        }


    };

} // Givaro

#include "givaro/givtruncdomain.inl"

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
