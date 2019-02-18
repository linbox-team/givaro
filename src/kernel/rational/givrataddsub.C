// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/rational/givrataddsub.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: M. Samama
// $Id: givrataddsub.C,v 1.4 2009-09-17 14:28:22 jgdumas Exp $
// ==========================================================================
// Description:
#include "givaro/givrational.h"

namespace Givaro {

    // ----------------------------------------- Rational::operator +
    Rational Rational::operator + (const Rational& r)  const
    {
        if (isZero(r)) return *this ;
        if (isZero(*this)) return r ;
        if (isInteger(*this) && isInteger(r))
            return Rational(num+r.num) ;

        if (Rational::flags == Rational::NoReduce)
            return Rational( num*r.den + r.num*den, den*r.den, 0) ;

        Integer d1 = gcd(den, r.den);
        if (d1 == 1)
            return Rational(num * r.den + r.num * den, den * r.den, 0);
        Integer t = num * (r.den / d1) + r.num * (den / d1);
        Integer d2 = gcd(t, d1);
        return Rational( t / d2, (den / d1) * (r.den / d2), 0 );
    }



    Rational& Rational::operator += (const Rational& r)
    {
        if (isZero(r)) return *this ;
        if (isZero(*this)) {
            num = r.num;
            den = r.den;
            return *this;
        }
        if (isInteger(*this) && isInteger(r)) {
            num += r.num;
            return *this;
        }
        if (Rational::flags == Rational::NoReduce) {
            num *= r.den;
            num += r.num * den;
            den *= r.den;
            return *this;
        }

        Integer d1 = gcd(den, r.den);
        if (d1 == 1) {
            num *= r.den;
            num += r.num * den;
            den *= r.den;
            return *this;
        }

        num *= (r.den / d1);
        num += ( r.num * (den / d1) );
        Integer d2 = gcd(num, d1);

        num /= d2;

        den /= d1;
        den *= r.den;
        den /= d2;

        return *this;
    }



    // ----------------------------------------- Rational::operator -
    Rational Rational::operator - (const Rational& r)  const
    {
        if (isZero(r)) return *this ;
        if (isZero(*this)) return Rational(-r.num,r.den, 0) ;
        if (isInteger(*this) && isInteger(r))
            return Rational(num-r.num) ;

        if (Rational::flags == Rational::NoReduce)
            return Rational( num*r.den - r.num*den, den*r.den, 0) ;

        Integer d1 = gcd(den, r.den);
        if (d1 == 1)
            return Rational(num * r.den - r.num * den, den * r.den, 0);
        Integer t = num * (r.den / d1) - r.num * (den / d1);
        Integer d2 = gcd(t, d1);
        return Rational( t / d2, (den / d1) * (r.den / d2), 0 );
    }



    Rational& Rational::operator -= (const Rational& r)
    {
        if (isZero(r)) return *this ;
        if (isZero(*this)) {
            num = -r.num;
            den = r.den; // GV Jeu  2 aoû 2018 17:30:53 CEST, a "-" was also there
            return *this;
        }
        if (isInteger(*this) && isInteger(r)) {
            num -= r.num;
            return *this;
        }
        if (Rational::flags == Rational::NoReduce) {
            num *= r.den;
            num -= r.num * den;
            den *= r.den;
            return *this;
        }

        Integer d1 = gcd(den, r.den);
        if (d1 == 1) {
            num *= r.den;
            num -= r.num * den;
            den *= r.den;
            return *this;
        }

        num *= (r.den / d1);
        num -= ( r.num * (den / d1) );
        Integer d2 = gcd(num, d1);

        num /= d2;

        den /= d1;
        den *= r.den;
        den /= d2;

        return *this;
    }




    // ----------------------------------------- Rational::operator -
    Rational Rational::operator - ()  const
    { return Rational (-num, den, 0) ; }

} // namespace Givaro
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
