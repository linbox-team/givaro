// =============================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier, J.-G. Dumas
// =============================================================

/*! @file givinteger.h
 * @ingroup integers
 * @brief Integer Domain class specialization.
 */

#ifndef __GIVARO_integer_H
#define __GIVARO_integer_H

#include "givaro-config.h"
#include "gmp++/gmp++.h"
#include "givaro/givbasictype.h"
#include "givaro/givinit.h"
#include "givaro/giverror.h"
#include "givaro/givranditer.h"
#include "givaro/random-integer.h"
#include "givaro/zring.h"
#include <string>

namespace Givaro {

    //------------------------------------ Class IntegerDom
    //! Integer Domain, Specialization of ZRing
    template<>
    class ZRing<Integer> : public UnparametricZRing<Integer> {
    public:
        using Self_t = ZRing<Integer>;
        using Parent_t = UnparametricZRing<Integer>;
        typedef Integer Rep;
        typedef Rep Element;

        using Parent_t::Parent_t; // inherit constructors

        int operator==( const Self_t&) const
        {
            return 1;
        }
        int operator!=( const Self_t&) const
        {
            return 0;
        }


        template<class XXX> XXX& convert(XXX& x, const Rep& a) const
        {
            return Caster<XXX,Rep>(x,a);
        }

        // -- Specialize arithmetic operators
        Rep& mul( Rep& r, const Rep& a, const Rep& b ) const
        {
            return Integer::mul(r,a,b);
        }
        Rep& mulin( Rep& r, const Rep& b ) const
        {
            return Integer::mulin(r,b);
        }
        Rep& div( Rep& r, const Rep& a, const Rep& b ) const
        {
            return Integer::div(r,a,b);
        }
        Rep& divin( Rep& r, const Rep& b ) const
        {
            return Integer::divin(r,b);
        }
        Rep& mod( Rep& r, const Rep& a, const Rep& b ) const
        {
            return Integer::mod(r,a,b);
        }
        Rep& modin( Rep& r,const Rep& b ) const
        {
            return Integer::modin(r,b);
        }
        Rep& add( Rep& r, const Rep& a, const Rep& b ) const
        {
            return Integer::add(r,a,b);
        }
        Rep& addin( Rep& r, const Rep& b ) const
        {
            return Integer::addin(r,b);
        }
        Rep& sub( Rep& r, const Rep& a, const Rep& b ) const
        {
            return Integer::sub(r,a,b);
        }
        Rep& subin( Rep& r, const Rep& b) const
        {
            return Integer::subin(r,b);
        }
        Rep& divmod( Rep& q, Rep& r, const Rep& a, const Rep& b ) const
        {
            return Integer::divmod(q,r,a,b);
        }
        Rep& divexact( Rep& q, const Rep& a, const Rep& b ) const
        {
            return Integer::divexact(q,a,b);
        }

        Rep& axpy( Rep& r, const Rep& a, const Rep& b, const Rep& c ) const
        {
            return Integer::axpy(r,a,b,c);
        }
        Rep& maxpy( Rep& r, const Rep& a, const Rep& b, const Rep& c ) const
        {
            return Integer::maxpy(r,a,b,c);
        }
        Rep& maxpyin( Rep& r, const Rep& a, const Rep& b) const
        {
            return Integer::maxpyin(r,a,b);
        }
        Rep& axmy( Rep& r, const Rep& a, const Rep& b, const Rep& c ) const
        {
            return Integer::axmy(r,a,b,c);
        }
        Rep& axpyin( Rep& r, const Rep& a, const Rep& b ) const
        {
            return Integer::axpyin(r,a,b);
        }
        Rep& axmyin( Rep& r, const Rep& a, const Rep& b ) const
        {
            return Integer::axmyin(r,a,b);
        }

        // -- unary methods
        Rep& neg( Rep& r, const Rep& a ) const
        {
            return Integer::neg(r,a);
        }
        Rep& negin( Rep& r ) const
        {
            return Integer::negin(r);
        }

        Element& quo (Element& q, const Element& a, const Element& b) const {return Integer::floor(q, a, b);}
        Element& rem (Element& r, const Element& a, const Element& b) const {return Integer::mod(r,a,b);}
        Element& quoin (Element& a, const Element& b) const{return quo(a,a,b);}
        Element& remin (Element& a, const Element& b) const {return modin(a,b);}
        void quoRem (Element& q, Element& r, const Element& a, const Element& b) const{Integer::divmod(q,r,a,b);}

        Element& logtwo(Element& z, const Element& x) const {
            return z = Element(uint64_t(x.bitsize() - 1));
        }

        // -- extended gcd  q = gcd(a,b) = u*a+v*b;
        Rep& gcd( Rep& g, Rep& u, Rep& v, const Rep& a, const Rep& b ) const
        {
            return ::Givaro::gcd(g, u, v, a, b);
        }
        Rep& gcd( Rep& g, const Rep& a, const Rep& b ) const
        {
            return ::Givaro::gcd(g, a, b);
        }
        Rep& gcdin( Rep& g, const Rep& a) const
        {
            Rep tmp(g);
            return ::Givaro::gcd(g, tmp, a);
        }
        Rep& lcm( Rep& l, const Rep& a, const Rep& b ) const
        {
            return ::Givaro::lcm(l, a, b);
        }
        Rep& lcmin( Rep& l, const Rep& a) const
        { Rep tmp(l); return lcm(l, tmp, a);
        }

        Element &dxgcd(Element &g, Element &s, Element &t, Element &u, Element &v, const Element &a, const Element &b) const
        {
            gcd(g,s,t,a,b);
            div(u,a,g);
            div(v,b,g);
            return g;
        }

        Rep& inv(Rep& u, const Rep& a, const Rep& b) const
        {
            return ::Givaro::inv(u,a,b);
        }
        Rep& invin(Rep& u, const Rep& b) const
        {
            return ::Givaro::invin(u,b);
        }

        Rep& invmod(Rep& u, const Rep& a, const Rep& b) const
        {
            return inv(u,a,b);
        }
        Rep& invmodin(Rep& u, const Rep& b) const
        {
            return invin(u,b);
        }

        bool ratrecon(Rep& num, Rep& den, const Rep& f, const Rep& m,
                      const Rep& numbound,
                      bool forcereduce = true, bool recurs = true) const ;
        bool RationalReconstruction(Rep&, Rep&, const Rep&, const Rep&) const;
        bool RationalReconstruction(Rep&, Rep&,
                                    const Rep&, const Rep&, const Rep&,
                                    bool = true, bool = true) const ;
        bool RationalReconstruction(Rep&, Rep&, const Rep&, const Rep&,
                                    const Rep&, const Rep&) const;
        // - return n^l
        Rep& pow(Rep& r, const Rep& n, const int64_t l) const
        {
            return r = ::Givaro::pow(n, l);
        }
        Rep& pow(Rep& r, const Rep& n, const uint64_t l) const
        {
            return r = ::Givaro::pow(n, l);
        }
        Rep& pow(Rep& r, const Rep& n, const int32_t l) const
        {
            return r = ::Givaro::pow(n, (int64_t)l);
        }
        Rep& pow(Rep& r, const Rep& n, const uint32_t l) const
        {
            return r = ::Givaro::pow(n, (uint64_t)l);
        }

        // - return square root of n
        Rep& sqrt(Rep& s, const Rep& n) const
        {
            return ::Givaro::sqrt(s,n);
        }
        Rep& sqrt(Rep& s, Rep& r, const Rep& n) const
        {
            return ::Givaro::sqrtrem(s,n, r);
        }
        // - base p logarithm of a
        int64_t logp(const Rep& a, const Rep& p) const
        {
            return ::Givaro::logp(a,p);
        }

        // - return n^e % m
        Rep& powmod(Rep& r, const Rep& n, const int64_t e, const Rep& m) const
        {
            return r = ::Givaro::powmod(n, e, m);
        }
        Rep& powmod(Rep& r, const Rep& n, const Rep& e, const Rep& m) const
        {
            return r = ::Givaro::powmod(n, e, m);
        }

        // - Misc
        uint64_t length (const Rep& a) const
        {
            return ::Givaro::length(a);
        }
        int sign   (const Rep& a) const
        {
            return ::Givaro::sign(a);
        }
        bool isZero (const Rep& a) const
        {
            return ::Givaro::isZero(a);
        }
        bool isOne  (const Rep& a) const
        {
            return ::Givaro::isOne(a);
        }
        bool isMOne  (const Rep& a) const
        {
            return ::Givaro::isMOne(a);
        }
        /// isUnit
        inline  bool isUnit (const Rep& x) const
        {
            return ::Givaro::isOne(x) || ::Givaro::isMOne(x);
        }

        /** @brief isDivisor (a, b)
         *  Test if b | a.
         */
        inline  bool isDivisor (const Element& a, const Element& b) const
        {
            Element r;
            if (::Givaro::isZero(b)) return ::Givaro::isZero(a);
            return ::Givaro::isZero(mod(r,a,b));
        }

        Element& abs(Element& x, const Element& a) const {
            return x=::Givaro::abs(a);
        }

        Element abs(const Element& a) const {
            return ::Givaro::abs(a);
        }

        int32_t compare(const Rep& a, const Rep& b) const {
            return ::Givaro::compare(a,b);
        }

        bool areEqual (const Rep& a, const Rep& b) const
        {
            return ::Givaro::compare(a,b) ==0;
        }
        bool areNEqual(const Rep& a, const Rep& b) const
        {
            return ::Givaro::compare(a,b) !=0;
        }
        bool areAssociates(const Element &x, const Element &y) const
        {
            return ::Givaro::abs(x) == ::Givaro::abs(y);
        }
        bool isgeq(const Rep& a, const Rep& b) const
        {
            return ::Givaro::compare(a,b) >= 0;
        }
        bool isleq(const Rep& a, const Rep& b) const
        {
            return ::Givaro::compare(a,b) <= 0;
        }
        bool isgeq(const int64_t b,const Rep& a ) const
        {
            return this->isgeq(Rep(b),a);
        }
        bool isleq(const int64_t b,const Rep& a ) const
        {
            return this->isleq(Rep(b),a);
        }
        bool isgeq(const Rep& a, const int64_t b) const
        {
            return this->isgeq(a,Rep(b));
        }
        bool isleq(const Rep& a, const int64_t b) const
        {
            return this->isleq(a,Rep(b));
        }
        bool isgt(const Rep& a, const Rep& b) const
        {
            return ::Givaro::compare(a,b) > 0;
        }
        bool islt(const Rep& a, const Rep& b) const
        {
            return ::Givaro::compare(a,b) < 0;
        }
        bool isgt(const int64_t b,const Rep& a ) const
        {
            return this->isgt(Rep(b),a);
        }
        bool islt(const int64_t b,const Rep& a ) const
        {
            return this->islt(Rep(b),a);
        }
        bool isgt(const Rep& a, const int64_t b) const
        {
            return this->isgt(a,Rep(b));
        }
        bool islt(const Rep& a, const int64_t b) const
        {
            return this->islt(a,Rep(b));
        }


#ifdef __GMP_PLUSPLUS__
        void seeding(unsigned long s = 0) const
        { Integer::seeding(s) ;
        }
#endif

        typedef RandomIntegerIterator<false,false> RandIter;
        typedef GeneralRingNonZeroRandIter<Self_t> NonZeroRandIter;


        template< class MyRandIter > Rep& random(MyRandIter&, Rep& r, long s = 1) const
        {
            return Integer::random(r,s);
        }
        template< class MyRandIter > Rep& random(MyRandIter&, Rep& r, const Rep& b) const
        {
            return Integer::random(r,b);
        }
        template< class MyRandIter > Rep& nonzerorandom(MyRandIter&, Rep& r, long s = 1) const
        {
            return Integer::nonzerorandom(r,s);
        }
        template< class MyRandIter > Rep& nonzerorandom (MyRandIter&,Rep& r, const Rep& b) const
        {
            return Integer::nonzerorandom(r,b);
        }

        // -- IO
        std::istream& read ( std::istream& i )
        {
            char ch;
            i >> std::ws >> ch;
            // JGD 22.03.03
            //    if (ch != 'I')
            //      GivError::throw_error(GivBadFormat("Self_t::read: bad signature domain"));
            return i;
        }
        std::ostream& write( std::ostream& o ) const
        {
            return o << "ZRing<Integer>";
        }

        using Parent_t::read;
        using Parent_t::write;
    };

    template<> struct DomainRandIter<ZRing<Integer>> {
        typedef ZRing<Integer>::RandIter RandIter;
    };

    typedef ZRing<Integer> IntegerDom;
    using IntegerDomain = ZRing<Integer>;


//     void RationalReconstruction(Integer& a, Integer& b,
//                                 const Integer& f, const Integer& m,
//                                 const Integer& k,
//                                 bool forcereduce, bool recursive );
//     bool ratrecon(Integer& num, Integer& den,
//                   const Integer& f, const Integer& m,
//                   const Integer& k,
//                   bool forcereduce = true, bool recurs = true );

} // Givaro

#endif //__GIVARO_integer_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
