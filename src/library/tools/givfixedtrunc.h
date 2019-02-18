// ===============================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <01 Apr 11 11:38:28 Jean-Guillaume.Dumas@imag.fr>
// Author: J-G. Dumas
// Description: Polynomials modulo X^{k+1}
// ===============================================================
#ifndef __GIVARO_fixed_trunc_domain_H
#include <givaro/givtruncdomain.h>

namespace Givaro {
    template <class Domain>
    class FixedTruncDom : public TruncDom<Domain> {
    public :
        // -- Self_t
        typedef          FixedTruncDom<Domain>	Self_t;
        // -- Father_t
        typedef          TruncDom<Domain>  	Father_t;
        typedef typename TruncDom<Domain>::Father_t PolDom;
        // -- Exported types
        typedef          Domain		       	Domain_t;
        typedef typename Domain::Element      Type_t;
        typedef typename Father_t::Polynomial_t	Polynomial_t;
        typedef typename Father_t::Storage_t	Storage_t;
        typedef          Storage_t                Rep;
        typedef          Storage_t                Element;

        Degree _deg;

        FixedTruncDom (const Domain& d, const Degree deg, const Indeter& X = Indeter() ) : Father_t(d,X), _deg(deg) {}
        FixedTruncDom (const Self_t& t) : Father_t(static_cast<const Father_t&>(t)),_deg(t._deg) {}
        FixedTruncDom (const Father_t& t, const Degree deg) : Father_t(t), _deg(deg) {}


        Degree& getModulus(Degree& d) const
        {
            return d=_deg;
        }

        Rep& init(Rep& p) const
        {
            return Father_t::init(p);
        }

        template<class XXX>
        Rep& init(Rep& p, const XXX &cste ) const
        {
            return Father_t::truncin(Father_t::init(p,cste),0,_deg);
        }

        // -- For polynomial = lcoeff X^deg
        template<class XXX>
        Rep& init (Rep& p, const Degree deg , const XXX& lcoeff) const
        {
            return Father_t::truncin(Father_t::init(p,deg,lcoeff),0,_deg);
        }
        //    F.assign(P[deg], lcoeff);
        Rep& assign (Rep& p, const Degree deg , const Type_t& lcoeff) const
        {
            return Father_t::truncin(Father_t::assign(p,deg,lcoeff),0,_deg);
        }
        // -- Assignment p = q
        Rep& assign( Rep& p, const Rep& q) const
        {
            return Father_t::assign(p,q);
        }

        Rep& assign(Rep& p, const Polynomial_t& r ) const
        {
            return Father_t::truncin(Father_t::assign(p,r),0,_deg);
        }

        Rep& assign(Rep& p, const Polynomial_t& r, const Degree v, const Degree d) const;


        Rep& mulin(Rep& p, const Degree& s) const
        {
            return Father_t::truncin(Father_t::mulin(p,s),0,_deg);
        }

        Rep& shiftin(Rep& p, const Degree& s) const
        {
            return Father_t::truncin(Father_t::mulin(p,s),0,_deg);
        }
        Rep& truncin(Rep& p, const Degree& v, const Degree& d) const
        {
            return Father_t::truncin(p,v,(d>_deg?_deg:d));
        }
        Rep& trunc(Rep& p, const Rep& R, const Degree& v, const Degree& d) const
        {
            return Father_t::trunc(p,R,v,(d>_deg?_deg:d));
        }


        // -- Arithmetics operators
        Rep& addin ( Rep& R, const Rep& P) const
        {
            return Father_t::addin(R, P, 0, _deg);
        }

        Rep& add ( Rep& res, const Rep& u, const Rep& v ) const
        {
            return Father_t::add(res, u, v, 0, _deg);
        }

        Rep& addin ( Rep& R, const Rep& P, const Degree& v, const Degree& d) const;

        Rep& add ( Rep& res, const Rep& u, const Rep& v, const Degree& val, const Degree& deg) const;

        Rep& sub ( Rep& res, const Rep& u, const Rep& v ) const
        {
            return Father_t::sub(res,u,v,0,_deg);
        }
        Rep& subin ( Rep& R, const Rep& P) const
        {
            return Father_t::subin(R,P,0,_deg);
        }
        Rep& sub ( Rep& R, const Rep& P, const Rep& Q, const Degree& v, const Degree& d) const;
        Rep& subin ( Rep& R, const Rep& P, const Degree& v, const Degree& d) const;

        Rep& mul ( Rep& res, const Rep& u, const Rep& v ) const
        {
            return Father_t::mul(res,u,v, 0, _deg);
        }
        Rep& mulin ( Rep& P, const Rep& Q ) const
        {
            return Father_t::mulin(P,Q, 0, _deg);
        }


        Rep& mul( Rep& r, const Rep& u, const Rep& v, const Degree& val, const Degree& deg) const;
        Rep& mulin( Rep& r, const Rep& v, const Degree& val, const Degree& deg) const;


        Rep& axpy  (Rep& r, const Rep& a, const Rep& x, const Rep& y) const
        {
            return Father_t::axpy(r,a,x,y, 0, _deg);
        }
        Rep& axpyin(Rep& r, const Rep& a, const Rep& x) const
        {
            return Father_t::axpyin(r,a,x, 0, _deg);
        }
        Rep& axpy  (Rep& r, const Rep& a, const Rep& x, const Rep& y, const Degree& val, const Degree& deg) const;

        Rep& axpyin  (Rep& r, const Rep& a, const Rep& x, const Degree& val, const Degree& deg) const;

        Rep& axmy  (Rep& r, const Rep& a, const Rep& x, const Rep& y) const
        {
            return Father_t::axmy(r,a,x,y, 0, _deg);
        }
        Rep& axmyin(Rep& r, const Rep& a, const Rep& x) const
        {
            return Father_t::axmyin(r,a,x, 0, _deg);
        }
        Rep& axmy  (Rep& r, const Rep& a, const Rep& x, const Rep& y, const Degree& val, const Degree& deg) const ;

        Rep& axmyin  (Rep& r, const Rep& a, const Rep& x, const Degree& val, const Degree& deg) const;


        Rep& maxpy  (Rep& r, const Rep& a, const Rep& x, const Rep& y) const
        {
            return Father_t::maxpy(r,a,x,y, 0, _deg);
        }
        Rep& maxpyin(Rep& r, const Rep& a, const Rep& x) const
        {
            return Father_t::maxpyin(r,a,x, 0, _deg);
        }
        Rep& maxpy  (Rep& r, const Rep& a, const Rep& x, const Rep& y, const Degree& val, const Degree& deg) const;

        Rep& maxpyin  (Rep& r, const Rep& a, const Rep& x, const Degree& val, const Degree& deg) const;



        Rep& invin ( Rep& q) const
        {
            Polynomial_t Xk; PolDom::init(Xk,_deg+1);
            Polynomial_t pq; Father_t::convert(pq,q);
            Polynomial_t t; PolDom::invmod(t,pq,Xk);
            return Father_t::assign(q,t);
        }
        Rep& inv( Rep& r, const Rep& u) const
        {
            return this->invin(Father_t::assign(r,u));
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




        // -- Random dense polynomial of degree d
        template< class RandIter > Rep& random(RandIter& g, Rep& r, Degree s) const
        {
            return Father_t::truncin(Father_t::random(g,r,s),0,_deg);
        }

        Rep& random(GivRandom& g, Rep& r, Degree s) const
        {
            return Father_t::truncin(Father_t::random(g,r,s),0,_deg);
        }

        Type_t& getEntry(Type_t& c, const Degree& i, const Rep& P) const
        {
            if (i>=_deg)
                return this->_domain.assign(c, this->_domain.zero);
            else
                return PolDom::getEntry(c,i,P.first);
        }

        Type_t& leadcoef (Type_t& c, const Rep& P) const
        {
            return PolDom::leadcoef(c,P.first);
        }

        std::ostream& write( std::ostream& o) const
        {
            return PolDom::write(o) << " mod " << this->_x << '^' << _deg;
        }

        std::istream& read ( std::istream& i, Rep& n) const
        {
            n.second=0;
            return Father_t::read(i,n.first);
        }
        std::ostream& write( std::ostream& o, const Rep& n) const
        {
            return Father_t::write(o,n);
        }

    };

} // Givaro

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
