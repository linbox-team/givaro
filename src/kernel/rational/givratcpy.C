// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/rational/givratcpy.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: M. Samama
// $Id: givratcpy.C,v 1.2 2009-09-17 14:28:23 jgdumas Exp $
// ==========================================================================
// Description:

#include "givaro/givrational.h"

namespace Givaro {

    Rational& Rational::logcpy (const Rational &r)
    {
        if (this == &r) return *this ;
        num.logcpy(r.num) ; den.logcpy(r.den) ;
        return *this ;
    }

    // same that Rational::logcpy function
    Rational& Rational::operator= (const Rational &r)
    {
        if (this == &r) return *this ;
        num.logcpy(r.num) ; den.logcpy(r.den) ;
        return *this ;
    }

    Rational& Rational::copy (const Rational &r)
    {
        if (this == &r) return *this ;
        num.copy(r.num) ; den.copy(r.den) ;
        return *this ;
    }

} // namespace Givaro

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
