// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/rational/givratmisc.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: M. Samama
// $Id: givratmisc.C,v 1.5 2009-10-01 09:07:36 jgdumas Exp $
// ==========================================================================
// Description:

#include "givaro/givrational.h"

namespace Givaro {

    //
    Rational& Rational::reduce()
    {
        Integer t = gcd(num, den);
        if (!isOne(t) )
        {
            num /= t;
            den /= t;
        }
        return *this;
    }

    const Integer trunc(const Rational &r)
    {
        return r.num / r.den;
    }

    const Integer floor(const Rational &x){return Integer::floor (x.num, x.den);}

    const Integer ceil(const Rational &x) {return Integer::ceil (x.num, x.den);}

    const Integer round(const Rational& x)  // GV Jeu  2 aoû 2018 16:38:58 CEST
    {
        Integer q;
        Integer r;
        Integer::divmod(q, r, abs(x.num), x.den);
        if (r!=0 && absCompare(r<<1, x.den) >= 0) // GV was < 0, and changed with abs
            q += 1;
        return (x.num<0) ? -q : q;

    }

    const Rational pow (const Rational& x, const int64_t y)
    {
        Rational r;
        if (y >= 0)
        {
            r.num = pow(x.num, (int64_t) y);
            r.den = pow(x.den, (int64_t) y);
        }
        else
        {
            r.den = pow(x.num, (int64_t) -y);
            r.num = pow(x.den, (int64_t) -y);
            if (sign(r.den) < 0)
            {
                r.num = -r.num ;
                r.den = -r.den ;
            }
        }
        return r;
    }

} // namespace Givaro
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
