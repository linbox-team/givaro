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

        /** @name Comparison operators.
         * @brief Compare with operators.
         */
        ///@{
        /** greater or equal.
         * @param l integer to be compared to
         */
        giv_all_inlined int32_t operator >= (const Integer & l) const;
        /** @overload Rational::operator>=(Rational) */
        giv_all_inlined int32_t operator >= (const int32_t l) const;
        /** @overload Rational::operator>=(Rational) */
        giv_all_inlined int32_t operator >= (const int64_t l) const;
        /** @overload Rational::operator>=(Rational) */
        giv_all_inlined int32_t operator >= (const uint64_t l) const;
        /** @overload Rational::operator>=(Rational) */
        giv_all_inlined int32_t operator >= (const uint32_t l) const;
        /** @overload Rational::operator>=(Rational) */
        giv_all_inlined int32_t operator >= (const double l) const;
        /** @overload Rational::operator>=(Rational) */
        giv_all_inlined int32_t operator >= (const float l) const;

        //! greater or equal.
        /// @param l,n rationals to compare
        giv_all_inlined friend int32_t operator >= (uint32_t l, const Rational& n);
        /** @overload Rational::operator>=(unsigned, Rational) */
        giv_all_inlined friend int32_t operator >= (float l, const Rational& n);
        /** @overload Rational::operator>=(unsigned, Rational) */
        giv_all_inlined friend int32_t operator >= (double l, const Rational& n);
        /** @overload Rational::operator>=(unsigned, Rational) */
        giv_all_inlined friend int32_t operator >= (int32_t l, const Rational& n);
        /** @overload Rational::operator>=(unsigned, Rational) */
        giv_all_inlined friend int32_t operator >= (int64_t l, const Rational& n);
        /** @overload Rational::operator>=(unsigned, Rational) */
        giv_all_inlined friend int32_t operator >= (uint64_t l, const Rational& n);


        //! less or equal
        /// @param l rational to be compared to
        giv_all_inlined int32_t operator <= (const Rational & l) const;
        /** @overload Rational::operator<=(Rational) */
        giv_all_inlined int32_t operator <= (const int32_t l) const;
        /** @overload Rational::operator<=(Rational) */
        giv_all_inlined int32_t operator <= (const int64_t l) const;
        /** @overload Rational::operator<=(Rational) */
        giv_all_inlined int32_t operator <= (const uint64_t l) const;
        /** @overload Rational::operator<=(Rational) */
        giv_all_inlined int32_t operator <= (const uint32_t l) const;
        /** @overload Rational::operator<=(Rational) */
        giv_all_inlined int32_t operator <= (const double l) const;
        /** @overload Rational::operator<=(Rational) */
        giv_all_inlined int32_t operator <= (const float l) const;

        //! less or equal
        /// @param l,n rationals to compare
        giv_all_inlined friend int32_t operator <= (uint32_t l, const Rational& n);
        /** @overload Rational::operator>=(unsigned, Rational) */
        giv_all_inlined friend int32_t operator <= (float l, const Rational& n);
        /** @overload Rational::operator>=(unsigned, Rational) */
        giv_all_inlined friend int32_t operator <= (double l, const Rational& n);
        /** @overload Rational::operator>=(unsigned, Rational) */
        giv_all_inlined friend int32_t operator <= (int32_t l, const Rational& n);
        /** @overload Rational::operator>=(unsigned, Rational) */
        giv_all_inlined friend int32_t operator <= (int64_t l, const Rational& n);
        /** @overload Rational::operator>=(unsigned, Rational) */
        giv_all_inlined friend int32_t operator <= (uint64_t l, const Rational& n);

        /*! operator != (not equal)
         * @param l rational
         * @return \c 1 iff l == this
         */
        giv_all_inlined int32_t operator != (const Integer & l) const;
        /** @overload Rational::operator!=(Rational) */
        giv_all_inlined int32_t operator != (const int32_t l) const;
        /** @overload Rational::operator!=(Rational) */
        giv_all_inlined int32_t operator != (const int64_t l) const;
        /** @overload Rational::operator!=(Rational) */
        giv_all_inlined int32_t operator != (const uint64_t l) const;
        /** @overload Rational::operator!=(Rational) */
        giv_all_inlined int32_t operator != (const uint32_t l) const;
        /** @overload Rational::operator!=(Rational) */
        giv_all_inlined int32_t operator != (const double l) const;
        /** @overload Rational::operator!=(Rational) */
        giv_all_inlined int32_t operator != (const float l) const;

        /*! operator != (not equal)
         * @param l,n rational
         * @return \c 1 iff l == n
         */
        giv_all_inlined friend int32_t operator != (uint32_t l, const Rational& n);
        /** @overload Rational::operator!=(unsigned, Rational) */
        giv_all_inlined friend int32_t operator != (float l, const Rational& n);
        /** @overload Rational::operator!=(unsigned, Rational) */
        giv_all_inlined friend int32_t operator != (double l, const Rational& n);
        /** @overload Rational::operator!=(unsigned, Rational) */
        giv_all_inlined friend int32_t operator != (int32_t l, const Rational& n);
        /** @overload Rational::operator!=(unsigned, Rational) */
        giv_all_inlined friend int32_t operator != (int64_t l, const Rational& n);
        /** @overload Rational::operator!=(unsigned, Rational) */
        giv_all_inlined friend int32_t operator != (uint64_t l, const Rational& n);


        //! Equality
        /// @param l rational to be compared to
        giv_all_inlined int32_t operator == (const Integer & l) const;
        /** @overload Rational::operator==(Rational) */
        giv_all_inlined int32_t operator == (const int32_t l) const;
        /** @overload Rational::operator==(Rational) */
        giv_all_inlined int32_t operator == (const int64_t l) const;
        /** @overload Rational::operator==(Rational) */
        giv_all_inlined int32_t operator == (const uint64_t l) const;
        /** @overload Rational::operator==(Rational) */
        giv_all_inlined int32_t operator == (const uint32_t l) const;
        /** @overload Rational::operator==(Rational) */
        giv_all_inlined int32_t operator == (const double l) const;
        /** @overload Rational::operator==(Rational) */
        giv_all_inlined int32_t operator == (const float l) const;

        //! Equality
        /// @param l,n rationals to compare
        giv_all_inlined friend int32_t operator == (uint32_t l, const Rational& n);
        /** @overload Rational::operator==(unsigned, Rational) */
        giv_all_inlined friend int32_t operator == (float l, const Rational& n);
        /** @overload Rational::operator==(unsigned, Rational) */
        giv_all_inlined friend int32_t operator == (double l, const Rational& n);
        /** @overload Rational::operator==(unsigned, Rational) */
        giv_all_inlined friend int32_t operator == (int32_t l, const Rational& n);
        /** @overload Rational::operator==(unsigned, Rational) */
        giv_all_inlined friend int32_t operator == (int64_t l, const Rational& n);
        /** @overload Rational::operator==(unsigned, Rational) */
        giv_all_inlined friend int32_t operator == (uint64_t l, const Rational& n);


        //! greater (strict)
        /// @param l rational to be compared to
        giv_all_inlined int32_t operator > (const Integer & l) const;
        /** @overload Rational::operator>(Rational) */
        giv_all_inlined int32_t operator > (const int32_t l) const;
        /** @overload Rational::operator>(Rational) */
        giv_all_inlined int32_t operator > (const int64_t l) const;
        /** @overload Rational::operator>(Rational) */
        giv_all_inlined int32_t operator > (const uint64_t l) const;
        /** @overload Rational::operator>(Rational) */
        giv_all_inlined int32_t operator > (const uint32_t l) const;
        /** @overload Rational::operator>(Rational) */
        giv_all_inlined int32_t operator > (const double l) const;
        /** @overload Rational::operator>(Rational) */
        giv_all_inlined int32_t operator > (const float l) const;

        //! greater (strict)
        /// @param l,n rationals to compare
        giv_all_inlined friend int32_t operator > (uint32_t l, const Rational& n);
        /** @overload Rational::operator>(unsigned, Rational) */
        giv_all_inlined friend int32_t operator > (float l, const Rational& n);
        /** @overload Rational::operator>(unsigned, Rational) */
        giv_all_inlined friend int32_t operator > (double l, const Rational& n);
        /** @overload Rational::operator>(unsigned, Rational) */
        giv_all_inlined friend int32_t operator > (int32_t l, const Rational& n);
        /** @overload Rational::operator>(unsigned, Rational) */
        giv_all_inlined friend int32_t operator > (int64_t l, const Rational& n);
        /** @overload Rational::operator>(unsigned, Rational) */
        giv_all_inlined friend int32_t operator > (uint64_t l, const Rational& n);

        //! less (strict)
        /// @param l rational to be compared to
        giv_all_inlined int32_t operator < (const Integer & l) const;
        /** @overload Rational::operator<(Rational) */
        giv_all_inlined int32_t operator < (const int32_t l) const;
        /** @overload Rational::operator<(Rational) */
        giv_all_inlined int32_t operator < (const int64_t l) const;
        /** @overload Rational::operator<(Rational) */
        giv_all_inlined int32_t operator < (const uint64_t l) const;
        /** @overload Rational::operator<(Rational) */
        giv_all_inlined int32_t operator < (const uint32_t l) const;
        /** @overload Rational::operator<(Rational) */
        giv_all_inlined int32_t operator < (const double l) const;
        /** @overload Rational::operator<(Rational) */
        giv_all_inlined int32_t operator < (const float l) const;

        //! less (strict)
        /// @param l,n rationals to compare
        giv_all_inlined friend int32_t operator < (uint32_t l, const Rational& n);
        /** @overload Rational::operator<(unsigned, Rational) */
        giv_all_inlined friend int32_t operator < (float l, const Rational& n);
        /** @overload Rational::operator<(unsigned, Rational) */
        giv_all_inlined friend int32_t operator < (double l, const Rational& n);
        /** @overload Rational::operator<(unsigned, Rational) */
        giv_all_inlined friend int32_t operator < (int32_t l, const Rational& n);
        /** @overload Rational::operator<(unsigned, Rational) */
        giv_all_inlined friend int32_t operator < (int64_t l, const Rational& n);
        /** @overload Rational::operator<(unsigned, Rational) */
        giv_all_inlined friend int32_t operator < (uint64_t l, const Rational& n);
        ///@}

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
        /** @brief num/den \equiv f modulo m, with a bound k on num
            forcereduce : must return a reduced fraction num/den
            recurs : tries to augment the given bound k if failure
        */
        bool ratrecon(const Integer& f, const Integer& m, const Integer& k,
                      bool forcereduce = true, bool recurs = false ) ;
    public:

        static bool ratrecon(
            Integer& num, Integer& den, const Integer& f, const Integer& m,
            const Integer& k, bool forcereduce = true, bool recurs = true );

        static bool RationalReconstruction(
            Integer& num, Integer& den, const Integer& f, const Integer& m);

        static bool RationalReconstruction(
            Integer& num, Integer& den, const Integer& f, const Integer& m,
            const Integer& numbound,
            bool forcereduce = true, bool recursive = true );

        // - exportation of the module
        static bool RationalReconstruction(
            Integer& num, Integer& den, const Integer& f, const Integer& m,
            const Integer& numbound, const Integer& denbound );

        static GivModule Module;
        // -- Cstor for Zero and One to delay initialization after the main
        Rational( givNoInit );
    }; // ----------------------------------- End of Class Rationalional

    extern std::istream& operator>> (std::istream& in, Rational& r) ;
} //namespace Givaro


#include "givaro/givrational.inl"

#endif // __GIVARO_rational_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
