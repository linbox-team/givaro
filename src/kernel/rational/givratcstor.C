// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/rational/givratcstor.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: M. Samama
// $Id: givratcstor.C,v 1.9 2011-02-04 14:50:07 briceboyer Exp $
// ==========================================================================
// Description:

#include <givaro/givconfig.h>
#include "givaro/givrational.h"
#include "givaro/giverror.h"
#include "givaro/givpower.h"
#include "givaro/udl.h"
#include <cmath>
#include <float.h>
#include <iostream>

#if !defined(__MWERKS__) && (defined(__GNUC__) && (__GNUC__ == 2))
#define __GIVARO_OLD_SSTREAM__
#include <strstream>
#else
// -- new interface for string stream
#include <sstream>
#endif

namespace Givaro {

    Rational::ReduceFlag Rational::flags = Rational::Reduce ;
    void Rational::SetReduce() { Rational::flags = Rational::Reduce ; }
    void Rational::SetNoReduce() { Rational::flags = Rational::NoReduce ; }

    // Explicit instanciation
    template double power(double x, unsigned int p) ;



    struct ieee {
#if HAVE_BIG_ENDIAN
        uint64_t negative:1;
        uint64_t exponent:11;
        uint64_t mantissa:52;
#else
        uint64_t mantissa:52;
        uint64_t exponent:11;
        uint64_t negative:1;
#endif
    };

    Rational::Rational(double x) {
        union {
            uint64_t l;
            ieee u;
            double d;
        } t; // temp

        t.d = x;
        if (t.u.exponent == 0) {
            // Denormal numbers
            num = (x<0.?-t.u.mantissa:t.u.mantissa);
            den = 1;
            *this/=Rational(Integer(1)<<1074);
        } else {
            const int64_t shift = 1075-t.u.exponent;
            t.u.exponent = 1076;
            if (shift > 0) {
                Integer tt( t.u.mantissa+static_cast<uint64_t>(4503599627370496_ui64) );
                num = (x<0.?-tt:tt);
                den = Integer(1)<<shift;
            } else {
                Integer tt( t.u.mantissa+static_cast<uint64_t>(4503599627370496_ui64));
                tt <<=(-shift);
                num = (x<0.?-tt:tt);
                den = 1;
            }
        }
        if (Rational::flags == Rational::Reduce) reduce();
    }


    //   ------------------------------ Rational(Neutral n )
    Rational::Rational(Neutral n ) : den(Integer::one)
    {
        if (n == Neutral::zero) {
            num = Integer::zero;
        }
        else { // n = one
            num = Integer::one;
        }
    }

    //   ------------------------------ Rational(int n)
    Rational::Rational(int32_t n ) : num(n), den(Integer::one)
    { }

    Rational::Rational(uint32_t n ) : num(n), den(Integer::one)
    { }


    //   ------------------------------ Rational(long n)
    Rational::Rational(int64_t n ) : num(n), den(Integer::one)
    { }


    //   ------------------------------ Rational(unsigned long n)
    Rational::Rational(uint64_t n ) : num(n), den(Integer::one)
    { }

    //   ------------------------------ Rational(unsigned long n, unsigned long d )
    Rational::Rational(uint64_t n, uint64_t d )
    {
        if (d == 0)
        {
            throw GivMathDivZero("[Rational::Rational]: null denominator of the rational.") ;
        }

        if (n == 0)
        {
            num = Integer::zero;
            den = Integer::one;
        }
        else
        {
            num = Integer(n);
            den = Integer(d);
        }
        reduce();
    }

    Rational::Rational(uint32_t n, uint32_t d ) : Rational( uint64_t(n), uint64_t(d) ) {}

    //   ------------------------------ Rational(int64_t n, int64_t d )
    Rational::Rational(int64_t n, int64_t d )
    {
        if (d == 0)
        {
            throw GivMathDivZero("[Rational::Rational]: null denominator of the rational.") ;
        }

        if (n == 0)
        {
            num = Integer::zero;
            den = Integer::one;
        }
        if (d > 0)
        {
            num = Integer(n);
            den = Integer(d);
        }
        else
        {
            num = Integer(-n);
            den = Integer(-d);
        }
        reduce();
    }
    Rational::Rational(int32_t n, int32_t d ) : Rational( int64_t(n), int64_t(d) ) {}


    //   ------------------------------ Rational(const char* s )
    Rational::Rational(const char* s )
    {
#ifdef __GIVARO_OLD_SSTREAM__
        std::istrstream input (s) ;
#else
        std::istringstream input (s) ;
#endif
        Rational r ;
        input >> r ;
        operator= (r) ;
    }


    //   ------------------------------ Rational(const Integer &n)
    Rational::Rational(const Integer &n) : den(Integer::one)
    {
        if (isZero(n))
        {
            num = Integer::zero;
        }
        else
        {
            num = n ;
        }
    }

    // ------------------------------ Rational(const Integer &n, const Integer &d)
    // If red == 1 then the rational is reduced (gcd computation!)
    Rational::Rational(const Integer &n, const Integer &d, int red)
    {
        if (isZero(d))
        {
            throw GivMathDivZero( "[Rational::Rational]: null denominator of the rational.") ;
        }

        if (isZero(n))
        {
            num = Integer::zero;
            den = Integer::one;
        }
        if (sign(d) > 0)
        {
            num = n ;
            den = d ;
        }
        else
        {
            num = -n;
            den = -d;
        }
        if (red == 1) reduce();
    }

    //   ------------------------------ Rational(const Rational &n)
    Rational::Rational(const Rational &r) : num(r.num), den(r.den)
    { }

    Rational::Rational( givNoInit ) : num(Integer::zero), den(Integer::one)
    { }

    // ------ Initialization module:
    //GivModule Rational::Module (Rational::Init,
    //                            Rational::End,
    //                            InitAfter(Integer::Module),
    //                            "Rational package") ;
    void Rational::Init(int* , char***)
    {}

    void Rational::End()
    {}

    const Rational Rational::one (1) ;
    const Rational Rational::mOne(-1) ;
    const Rational Rational::zero(0) ;
} // namespace Givaro

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
