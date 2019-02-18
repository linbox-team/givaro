// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givindeter.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givindeter.h,v 1.6 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================
/** @file givindeter.h
 * @ingroup poly1
 * @brief indeterminates for polynomial manipulation
 */
#ifndef __GIVARO_indeter_H
#define __GIVARO_indeter_H

#include <iostream>
#include <string>

namespace Givaro {

    //! Indeterminate
    class Indeter {
    public :

        // -- Cstor: recopy the string
        Indeter(const std::string & x="") : name(x){}
        // -- Cstor: recopy the string
        Indeter(const char * x) : name(x){}
        // -- Cstor: recopy the string
        Indeter(const char c) : name(1U,c){}
        // -- Cstor of recopy
        Indeter(const Indeter& s): name(s.name) {}

        // -- Dstor
        ~Indeter(){}

        // -- assignement
        Indeter& operator= (const Indeter& s);

        // -- Comparizon operators:
        // all comparizons are based on this virtual method,
        // which returns : -1 iff *this < b, 0 iff *this == b and
        // +1 else. This comparizon method gives the natural order
        // for multivariate polynomials.
        int compare(const Indeter& b)  const;

        // -- methods
        friend std::ostream& operator<< (std::ostream& o, const Indeter& X);
        friend std::istream& operator>> (std::istream& o, Indeter& X);

    protected:
        std::string name;
    };


    //! @bug put elsewere. Inline members functions :
    inline int operator==(const Indeter& i1, const Indeter &i2)
    { return i1.compare(i2) ==0; }

    inline int operator!=(const Indeter& i1, const Indeter &i2)
    { return i1.compare(i2) !=0; }

    inline int operator<= (const Indeter& i1, const Indeter &i2)
    { return i1.compare(i2) <=0; }

    inline int operator<  (const Indeter& i1, const Indeter &i2)
    { return i1.compare(i2) <0; }

    inline int operator>= (const Indeter& i1, const Indeter &i2)
    { return i1.compare(i2) >=0; }

    inline int operator>  (const Indeter& i1, const Indeter &i2)
    { return i1.compare(i2) >0; }

} // Givaro
#endif // __GIVARO_indeter_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
