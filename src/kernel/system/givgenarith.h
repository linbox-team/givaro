// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/system/givgenarith.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givgenarith.h,v 1.5 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================
/** @file givgenarith.h
 * @ingroup system
 * @brief Domain definition for basic type of the language.
 */
#ifndef __GIVARO_genarith_H
#define __GIVARO_genarith_H

#include "givaro/givbasictype.h"
namespace Givaro {
    //! give a name for /read/write
    template<class T>
    struct __givdom_trait_name { enum { val = '?' }; };

#define __GIVARO_SPEC_NAME( type, name ) \
    struct __givdom_trait_name<type> { enum { val = name }; } ;

    //! Base Domain
    template<class T>
    class BaseDomain {
    public:
        typedef T	Rep;
        enum { size_rep = sizeof(Rep) };

        const Rep zero;
        const Rep one;
        const Rep mone;

        BaseDomain() : zero((T)0), one((T)1), mone((T)-1) {}

        // -- Assignment of domain: default operator=
        // -- Comparizon of domain:
        int operator==( const BaseDomain<T>& BC) const { return 1;}
        int operator!=( const BaseDomain<T>& BC) const { return 0;}

        // -- Misc
        void init  ( Rep& r ) const {};
        void init  ( Rep& r, const Rep a ) const { r = a; };
        void assign( Rep& r, const Rep a ) const { r = a; };

        // -- Comparizon of Rep
        int isZero( const Rep r ) const { return r ==zero; };
        int isOne ( const Rep r ) const { return r ==one; };
        int isMOne ( const Rep r ) const { return r ==mone; };
        int areEqual ( const Rep r, const Rep a ) const { r == a; };
        int areNEqual( const Rep r, const Rep a ) const { r != a; };

        // -- Arithmetic
        void mul( Rep& r, const Rep a, const Rep b ) const { r = a * b; };
        void div( Rep& r, const Rep a, const Rep b ) const { r = a / b; };
        void mod( Rep& r, const Rep a, const Rep b ) const { r = a % b; };
        void add( Rep& r, const Rep a, const Rep b ) const { r = a + b; };
        void sub( Rep& r, const Rep a, const Rep b ) const { r = a - b; };

        void mulin( Rep& r, const Rep a) const { r *= a; };
        void divin( Rep& r, const Rep a) const { r /= a; };
        void modin( Rep& r, const Rep a) const { r %= a; };
        void addin( Rep& r, const Rep a) const { r += a; };
        void subin( Rep& r, const Rep a) const { r -= a; };

        void axpy( Rep& r, const Rep a, const Rep b, const Rep c ) const
        { r = a * b + c; };
        void axpyin( Rep& r, const Rep a, const Rep b ) const { r += a * b; };
        void axmy( Rep& r, const Rep a, const Rep b, const Rep c ) const
        { r = a * b - c; };
        void axmyin( Rep& r, const Rep a, const Rep b ) const { r = a * b - r; };

        // -- unary methods
        void neg( Rep& r, const Rep a ) const { r = -a; };
        void inv( Rep& r, const Rep a ) const { r = 1/a; };
        void negin( Rep& r ) const { r = -r; };
        void invin( Rep& r ) const { r = 1/r; };

        // -- IO methods: of the object domain
        ostream& write( ostream& s ) const
        { return s << (char)__givdom_trait_name<T>::val; }
        istream& read ( istream& s )
        {
            char c; s >> std::ws >> c;
            if (c != (char)__givdom_trait_name<T>::val)
                GivError::throw_error(GivBadFormat("Bad domain"));
        }

        // -- IO methods: of the element of the domain
        ostream& write( ostream& s, const Rep& r ) const
        { return s << r; }
        istream& read ( istream& s, Rep& r ) const
        { return s >> r; }
    };

    //! char dom
    typedef BaseDomain<char>    CharDom;
    __GIVARO_SPEC_NAME(char, 'c')
    //! short dom
    typedef BaseDomain<short>   ShortDom;
    __GIVARO_SPEC_NAME(short, 's')
    //! int dom
    typedef BaseDomain<int>     IntDom;
    __GIVARO_SPEC_NAME(int, 'i')
    //! long dom
    typedef BaseDomain<long>    LongDom;
    __GIVARO_SPEC_NAME(long, 'l')
    //! float dom
    typedef BaseDomain<float>   FloatDom;
    __GIVARO_SPEC_NAME(float, 'f')
    //! double dom
    typedef BaseDomain<double>  DoubleDom;
    __GIVARO_SPEC_NAME(double, 'd')

} // namespace Givaro

#endif // __GIVARO_genarith_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
