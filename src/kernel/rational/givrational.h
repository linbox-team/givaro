// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/rational/givrational.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: M. Samama, T. Gautier
// $Id: givrational.h,v 1.13 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================
/*! @file givrational.h
 * @ingroup rational
 * @brief Rationals (and domain),
 *        composed of an integer (numerator),
 *        and a positive integer (denominator)
 * NO DOC.
 */
#ifndef __GIVARO_rational_H
#define __GIVARO_rational_H
// #define __GIVARO_GMPplusplus_rational_H

#include "givaro/givinteger.h"
#include "givaro/givmodule.h"
#include "givaro/ring-interface.h"

namespace Givaro {

    // ----------------------------------- Functions Rational

    class Rational ;
    int compare(const Rational& a, const Rational& b) ;
    int absCompare(const Rational& a, const Rational& b) ;
    const Rational pow(const Rational &r, const int64_t l);
    const Integer floor (const Rational &r) ;
    const Integer ceil  (const Rational &r) ;
    const Integer round (const Rational &r) ;
    const Integer trunc (const Rational &r) ;
    const Rational abs  (const Rational &r) ;
    const Rational pow (const Rational& n, uint l);
    const Rational pow (const Rational& n, uint64_t l);
    uint64_t length (const Rational& r) ;
    int sign   (const Rational& r) ;
    int isZero (const Rational& r) ;
    int isOne  (const Rational& r) ;
    int isMOne  (const Rational& r) ;
    int isInteger(const Rational& r);

    template<class RatElement>
    class QField;


    //! Rationals. No doc.
    class Rational {

    public :
        // Cstor et dstor
        Rational(Neutral n = Neutral::zero) ;
        Rational(int32_t n) ;
        Rational(uint32_t n) ;
        Rational(int64_t n) ;
        Rational(uint64_t n) ;
        Rational(int64_t n, int64_t d ) ;
        Rational(uint64_t n, uint64_t d ) ;
        Rational(int32_t n, int32_t d ) ;
        Rational(uint32_t n, uint32_t d ) ;
        Rational(double x) ;
        Rational(const char* s) ;
        Rational(const Integer& n) ;
        Rational(const Integer& n, const Integer& d, int red = 1 ) ;
        // Rational number reconstruction
        /*! @brief Rational number reconstruction.
         * \f$ num/den \equiv f \mod m\f$, with \f$|num|<k\f$ and \f$0 < |den| \leq f/kf\f$
         * @bib
         * - von zur Gathen & Gerhard <i>Modern Computer Algebra</i>, 5.10, Cambridge Univ. Press 1999]
         */

        Rational(const Integer& f, const Integer& m, const Integer& k, bool recurs = false ) ;
        Rational(const Rational&) ;
        //~Rational();

        // Predefined cstes
        static const Rational zero ;
        static const Rational one ;
        static const Rational mOne ;

        // Logical and physical copies
        Rational& operator = (const Rational& );
        Rational& logcpy (const Rational& ) ;
        Rational& copy (const Rational& ) ;

        //------------------Equalities and inequalities between rationals
        friend int compare(const Rational& a, const Rational& b) ;
        friend int absCompare(const Rational& a, const Rational& b) ;


        //----------------Elementary arithmetic between Rational
        Rational operator + (const Rational& r) const ;
        Rational operator - (const Rational& r) const ;
        Rational operator - () const ;
        Rational operator + () const ;
        Rational operator * (const Rational& r) const ;
        Rational operator / (const Rational &r) const ;
        Rational& operator += (const Rational& r) ;
        Rational& operator -= (const Rational& r) ;
        Rational& operator *= (const Rational& r) ;
        Rational& operator /= (const Rational &r) ;

        Integer operator % (const Integer &r) const;

        //-----------------------------------------Arithmetic functions
        friend const Rational pow(const Rational &r, const int64_t l);

        //-----------------------------------------Miscellaneous
        friend const Integer floor (const Rational &r) ;
        friend const Integer ceil  (const Rational &r) ;
        friend const Integer round (const Rational &r) ;
        friend const Integer trunc (const Rational &r) ;

        inline friend const Rational abs  (const Rational &r) ;

        friend const Rational pow (const Rational& n, uint32_t l) {
            Rational r;
            r.num = ::Givaro::pow(n.num, l); r.den =  ::Givaro::pow(n.den, l);
            return r;
        }

        friend const Rational pow (const Rational& n, uint64_t l) {
            Rational r;
            r.num =  ::Givaro::pow(n.num, l); r.den =  ::Givaro::pow(n.den, l);
            return r;
        }

        const Integer nume() const ;
        const Integer deno() const ;
        inline friend uint64_t length (const Rational& r) ;
        inline friend int sign   (const Rational& r) ;
        inline friend int isZero (const Rational& r) ;
        inline friend int isOne  (const Rational& r) ;
        inline friend int isMOne  (const Rational& r) ;
        inline friend int isInteger(const Rational& r);

        std::ostream& print ( std::ostream& o ) const ;

        inline Rational reduce(const Rational& R) const ;

        static void SetReduce() ;
        static void SetNoReduce() ;

        // -- Cast operators
        operator short() const { return (short)(int) *this; }
        operator uint16_t() const { return (uint16_t) (uint32_t) *this; }
        operator uint8_t() const { return (uint8_t)(uint32_t) *this; }
        operator uint32_t() const { return (uint32_t) (this->num/this->den); }
        operator int() const  { return (int) (this->num/this->den); }
        operator signed char() const { return (signed char) (int) *this; }
        operator uint64_t() const { return (uint64_t) (this->num/this->den); }
        operator int64_t() const { return (int64_t) (this->num/this->den); }
        operator std::string() const { return std::string(this->num)+'/'+std::string(this->den); }
        operator float() const { return ((float)this->num)/((float)this->den); }
        operator double() const { return ((double)this->num)/((double)this->den); }

    protected: // Internal Representation : num/den
        Integer num, den;

    public:
        enum ReduceFlag { Reduce = 0x1, NoReduce = 0x0 } ;
    protected:
        static ReduceFlag flags ;    // flags that indicates is reduction is done or not
        // by default = Reduce
        Rational& reduce() ;
        // -- module initialization
        static void Init(int* argc, char***argv);
        static void End();
        friend class GivModule;
        friend class QField<Rational>;
        // Rational number reconstruction
        bool ratrecon(const Integer& f, const Integer& m, const Integer& k, bool recurs = false ) ;

    public:
        // - exportation of the module
        static GivModule Module;
        // -- Cstor for Zero and One to delay initialization after the main
        Rational( givNoInit );
    }; // ----------------------------------- End of Class Rationalional

    extern std::istream& operator>> (std::istream& in, Rational& r) ;
}


#include "givaro/givrational.inl"

namespace Givaro {

    //! Rational Domain
    template<>
    class QField<Rational> : public FieldInterface<Rational> {
    public:
        using Self_t = QField<Element>;
        typedef Rational Element;
        typedef Rational Rep;
        typedef uint64_t Residu_t; // type for cardinality()

        // -- Cstor
        QField() : one(1), mOne(-one), zero(0) {}
        template<class X> QField(const X& x) : one(1), mOne(-one),zero(0) {}

        int operator==( const QField<Element>& ) const { return 1;}
        int operator!=( const QField<Element>& ) const { return 0;}

        // -- Constants
        const Element one;
        const Element mOne;
        const Element zero;

        uint64_t characteristic() const { return 0U; }
        Integer& characteristic(Integer& p) const { return p=characteristic();}
        uint64_t cardinality() const { return 0U; }
        Integer& cardinality(Integer& p) const { return p=cardinality();}

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
        Rep& inv( Rep& r, const Rep& a ) const {
            const int snum( sign(a.num) );
#ifdef GIVARO_DEBUG
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
        Rep& negin( Rep& r ) const { Integer::negin(r.num); return r; }
        Rep& invin( Rep& r ) const {
            const int snum( sign(r.num) );
#ifdef GIVARO_DEBUG
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

#endif // __GIVARO_rational_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
