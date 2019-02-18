// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/rational/givratcompare.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: M. Samama
// $Id: givratcompare.C,v 1.3 2009-09-17 14:28:23 jgdumas Exp $
// ==========================================================================
// Description:
#include "givaro/givrational.h"

namespace Givaro {

    // -------------------------------------------- compare
    // returns 1 if a > b, 0 if a == b and -1 otherwise.
    int compare (const Rational& a, const Rational &b)
    {
        if (isZero(a.num) && isZero(b.num))
            return 0 ;
        if (isZero(a.num))
            return -sign(b.num);
        if (isZero(b.num))
            return sign(a.num);
        if (sign(a.num) != sign(b.num))
            return (sign(a.num) == -1 ? -1 : 1 )  ;
        if (sign(a.num) > 0)
            return absCompare(a,b) ;
        else return -absCompare(a,b) ;
    }

    // -------------------------------------------- absCompare
    int absCompare (const Rational& a, const Rational& b)
    {
        int cnum = absCompare(a.num, b.num) ;
        int cden = absCompare(a.den, b.den) ;

        if ( (cnum == -1) && (cden == 1) )
            return -1 ;
        if ( (cnum == 1) && (cden == -1) )
            return 1 ;
        if (cnum == 0)
            return -cden ;
        if (cden == 0)
            return cnum ;

        Integer p1 = a.num * b.den;
        Integer p2 = a.den * b.num;
        return absCompare(p1,p2);
    }

} // namespace Givaro
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
