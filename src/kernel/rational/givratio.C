// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/rational/givratio.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: M. Samama
// $Id: givratio.C,v 1.2 2009-09-17 14:28:23 jgdumas Exp $
// ==========================================================================
// Description:

#include <iostream>
#include "givaro/givrational.h"

namespace Givaro {

    std::ostream& Rational::print(std::ostream& s) const
    {
        if (den > 1) {
            s << num << "/" << den ;
        }
        else
            s << num;
        return s;
    }


    std::istream& operator>> (std::istream& in, Rational& r)
    {
        Integer num ;
        Integer den = 1;
        char ch ;

        in >> num ;
        if (!in.good() || in.eof()) {
            r = Rational(num) ;
            return in ;
        }

        if (in) {
            in.get(ch) ;
            if (in.eof()) {
                r = Rational(num) ;
                return in ;
            }
            while ((ch==' ') && (in)) in.get(ch) ;
            if (ch == '/') {
                // We get denominator
                in >> den ;
            }
            else in.putback(ch) ;
        } ;
        r = Rational(num,den) ;
        return in ;
    }

} // namespace Givaro
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
