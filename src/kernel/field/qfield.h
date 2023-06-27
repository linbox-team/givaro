// ==========================================================================
// $Source: qfield.h
// Copyright(c)'2019 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: M. Samama, T. Gautier
// $Id: givrational.h,v 1.13 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================
/*! @file qfield.h
 * @ingroup QField
 * @brief Specialization of Rational Domain
 * NO DOC.
 */
#ifndef __GIVARO_qfield_H
#define __GIVARO_qfield_H

#include "givaro/givrational.h"
#include "givaro/zring.h"

namespace Givaro {

    template<class RatElement>
    class QField;

    //! Rational Domain
    template<>
    class QField<Rational> : public FieldInterface<Rational> {
    public:
        using Self_t = QField<Element>;
        typedef Rational Element;
        typedef Rational Rep;
        typedef int64_t Residu_t; // type for cardinality()

        // -- Cstor
        QField() : one(1), mOne(-one), zero(0) {}
        template<class X> QField(const X& x) : one(1), mOne(-one),zero(0) {}

        int operator==( const QField<Element>& ) const { return 1;}
        int operator!=( const QField<Element>& ) const { return 0;}

        // -- Constants
        const Element one;
        const Element mOne;
        const Element zero;

        Residu_t characteristic() const { return 0; }
        Residu_t cardinality() const { return 0; }
        template<typename T> T& cardinality(T& c) const { return c = static_cast<T>(0); }
        template<typename T> T& characteristic(T& c) const { return c = static_cast<T>(0); }

        static inline Residu_t maxCardinality() { return -1; }
        static inline Residu_t minCardinality() { return 2; }

        // -- assignement
        Rep& init( Rep& a ) const{ return a; }
        Rep& init( Rep& a, const Integer& n, const Integer& d) const{ return a=Rational(n,d); }
        template<class XXX> Rep& init(Rep& r, const XXX& x) const {
            return Caster<Rep,XXX>(r,x);
        }
        template<class XXX> XXX& convert(XXX& x, const Rep& a) const {
            return Caster<XXX,Rep>(x,a);
        }

        Rep& assign( Rep& a, const Rep& b) const { return a = b ; }
        // -- integers operators
        Integer& get_num(Integer& n, const Element& r) const { return n=r.nume();}
        Integer& get_den(Integer& d, const Element& r) const { return d=r.deno();}

        // -- arithmetic operators
        Rep& mul( Rep& r, const Rep& a, const Rep& b ) const { return r = a * b; };
        Rep& div( Rep& r, const Rep& a, const Rep& b ) const { return r = a / b; };
        Rep& add( Rep& r, const Rep& a, const Rep& b ) const { return r = a + b; };
        Rep& sub( Rep& r, const Rep& a, const Rep& b ) const { return r = a - b; };

        Rep& mulin( Rep& r, const Rep& a) const { return r *= a; };
        Rep& divin( Rep& r, const Rep& a) const { return r /= a; };
        Rep& addin( Rep& r, const Rep& a) const { return r += a; };
        Rep& subin( Rep& r, const Rep& a) const { return r -= a; };

        Rep& axpy( Rep& r, const Rep& a, const Rep& b, const Rep& c ) const
        { return r = a * b + c; };
        Rep& axpyin( Rep& r, const Rep& a, const Rep& b ) const
        { return r += a * b; };
        Rep& maxpy( Rep& r, const Rep& a, const Rep& b, const Rep& c ) const
        { return r = c - a * b; };
        Rep& axmy( Rep& r, const Rep& a, const Rep& b, const Rep& c ) const
        { return r = a * b - c; };
        Rep& axmyin( Rep& r, const Rep& a, const Rep& b ) const
        { return r = a * b - r ; };
        Rep& maxpyin( Rep& r, const Rep& a, const Rep& b ) const
        { return r -= a * b; };

        // -- unary methods
        Rep& neg( Rep& r, const Rep& a ) const {
            Integer::neg(r.num,a.num);
            r.den=a.den; return r; }

        Rep& negin( Rep& r ) const { Integer::negin(r.num); return r; }

        Rep& inv( Rep& r, const Rep& a ) const {
            const int snum( sign(a.num) );
#ifdef __GIVARO_DEBUG
            if (snum == 0)
                throw GivMathDivZero("*** Error: division by zero, in operator Rational::inv in givrational.h") ;
#endif
            r.num=a.den; r.den=a.num;
            if (snum < 0) {
                Integer::negin(r.num);
                Integer::negin(r.den);
            }
            return r;
        }
        Rep& invin( Rep& r ) const {
            const int snum( sign(r.num) );
#ifdef __GIVARO_DEBUG
            if (snum == 0)
                throw GivMathDivZero("*** Error: division by zero, in operator Rational::invin in givrational.h") ;
#endif
            std::swap(r.num,r.den);
            if (snum < 0) {
                Integer::negin(r.num);
                Integer::negin(r.den);
            }
            return r;
        }

        // - return n^l
        Rep& pow(Rep& r, const Rep& n, const uint64_t l) const { return r =  ::Givaro::pow(n, l); }
        Rep& pow(Rep& r, const Rep& n, const uint32_t l) const { return r =  ::Givaro::pow(n, l); }


        // - Rational number reconstruction
        Rep& ratrecon(Rep& r, const Integer& f, const Integer& m, const Integer& k, bool recurs = false) const {
            return r = Rational(f,m,k,recurs);
        }
        Rep& ratrecon(Rep& r, const Integer& f, const Integer& m, bool recurs=true) const {
            return r = Rational(f,m, ::Givaro::sqrt(m),recurs);
        }


        // - Misc
        size_t length (const Rep& a) const { return  ::Givaro::length(a); }
        int sign    (const Rep& a) const { return  ::Givaro::sign(a); }
        bool isOne   (const Rep& a) const { return compare(a, one) ==0; }
        bool isMOne   (const Rep& a) const { return compare(a, mOne) ==0; }
        bool isUnit   (const Rep& a) const { return !isZero(a); }
        bool isZero  (const Rep& a) const { return compare(a, zero) ==0; }
        bool areEqual (const Rep& a, const Rep& b) const { return compare(a, b) ==0; }
        int areNEqual(const Rep& a, const Rep& b) const { return compare(a, b) !=0; }
        typedef GeneralRingRandIter<Self_t> RandIter;
        typedef GeneralRingNonZeroRandIter<Self_t> NonZeroRandIter;

        template< class MyRandIter > Rep& random(MyRandIter& g, Rep& r, int64_t s = 1) const { return r=Rational(Integer::random(s), Integer::nonzerorandom(s)); }
        template< class MyRandIter > Rep& random(MyRandIter& g, Rep& r, const Rep& b) const { Integer rnum,rden; Integer::random(rnum,b.nume()); Integer::nonzerorandom(rden,b.deno()); return r=Rational(rnum,rden); }
        template< class MyRandIter > Rep& nonzerorandom(MyRandIter& g, Rep& r, int64_t s = 1) const { return r=Rational(Integer::nonzerorandom(s), Integer::nonzerorandom(s)); }
        template< class MyRandIter > Rep& nonzerorandom (MyRandIter& g,Rep& r, const Rep& b) const { Integer rnum,rden; Integer::nonzerorandom(rnum,b.nume()); Integer::nonzerorandom(rden,b.deno()); return r=Rational(rnum,rden); }


        // -- IO
        // -- IO
        std::istream& read ( std::istream& i )
        { char ch;
            i >> std::ws >> ch;
            if (ch != 'R')
                GivError::throw_error(GivBadFormat("QField<Element>::read: bad signature domain"));
            return i;
        }
        std::ostream& write( std::ostream& o ) const { return o << 'R'; }
        std::istream& read ( std::istream& i, Rep& n) const { return i >> n; }
        std::ostream& write( std::ostream& o, const Rep& n) const { return n.print(o); }
    };

} //namespace Givaro

#endif // __GIVARO_qfield_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
