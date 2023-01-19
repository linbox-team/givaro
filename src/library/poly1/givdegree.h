// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givdegree.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givdegree.h,v 1.7 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================
/** @file givdegree.h
 * @ingroup poly1
 * @brief NO DOC
 * opaque class for Degree of polynomial. Degree of polynomial
 * 0 is Degree::deginfty with value DEGPOLYZERO.
 *
 */

#ifndef __GIVARO_poly1degree_H
#define __GIVARO_poly1degree_H

#include <cstdint>

#include <iostream>

namespace Givaro {
    //! Degree type for polynomials
    class Degree {
    public:
        typedef int64_t value_type;

        enum { DEGPOLYZERO =-1};
        Degree(int64_t a = DEGPOLYZERO): _deg(a<0?int64_t(DEGPOLYZERO):a) {}

        // -- Degree of zero polynomial
        static const int64_t deginfty;

        // -- cvrt
        int64_t value() const { return _deg; }

        // -- Basic arithmetic:
        Degree operator+( const Degree& d) const { return Degree(_deg+d._deg); }
        Degree operator-( const Degree& d) const { return Degree(_deg-d._deg); }
        Degree operator*( const uint64_t& e) const { return Degree(_deg*e); }
        Degree operator/( const uint64_t& e) const { return Degree(_deg/e); }
        Degree& operator+=( const Degree& d) { _deg+=d._deg; return *this; }
        Degree& operator-=( const Degree& d) { _deg-=d._deg; return *this; }
        Degree& operator*=( const uint64_t& e) { _deg*=e; return *this; }
        Degree& operator/=( const uint64_t& e) { _deg/=e; return *this; }
        Degree operator<<( const int i) const { return Degree(_deg<<i); }
        Degree operator>>( const int i) const { return Degree(_deg>>i); }
        Degree& operator <<=( const int i) { _deg<<=i; return *this;}
        Degree& operator >>=( const int i) { _deg>>=i; return *this;}
        int64_t operator++() { return ++_deg; }
        int64_t operator--() { return --_deg; }
        int64_t operator++(int) { return _deg++; }
        int64_t operator--(int) { return _deg--; }

        // -- Comparizon:
        int operator==( const Degree& d) const { return _deg == d._deg; }
        int operator!=( const Degree& d) const { return _deg != d._deg; }
        int operator<=( const Degree& d) const { return _deg <= d._deg; }
        int operator< ( const Degree& d) const { return _deg <  d._deg; }
        int operator>=( const Degree& d) const { return _deg >= d._deg; }
        int operator> ( const Degree& d) const { return _deg >  d._deg; }
        int operator==( const int64_t& d) const { return _deg == d; }
        int operator!=( const int64_t& d) const { return _deg != d; }
        int operator<=( const int64_t& d) const { return _deg <= d; }
        int operator< ( const int64_t& d) const { return _deg <  d; }
        int operator>=( const int64_t& d) const { return _deg >= d; }
        int operator> ( const int64_t& d) const { return _deg >  d; }

        // -- methods
        friend std::ostream& operator<< (std::ostream& o, const Degree& d) { return o << d._deg; }
        friend std::istream& operator>> (std::istream& i, Degree& d) { return i >> d._deg; }


    public:
        int64_t _deg;
    };

    //! value
    inline int64_t value(const Degree& d) { return d.value(); }
} // Givaro

#endif // __GIVARO_poly1degree_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
