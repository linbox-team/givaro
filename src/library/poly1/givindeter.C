// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givindeter.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givindeter.C,v 1.5 2009-09-17 14:28:23 jgdumas Exp $
// ==========================================================================
// Description:

#include <iostream>
#include <string.h>
#include "givaro/givindeter.h"

namespace Givaro {

    Indeter& Indeter::operator=( const Indeter& s )
    {
        name = s.name;
        return *this;
    }

    int Indeter::compare(const Indeter& b)  const
    {
        return name.compare(b.name);
    }

    std::ostream& operator<< (std::ostream& o, const Indeter& X)
    {
        //   return o << '[' << X.name.baseptr() << ']';
        return o << X.name ;
    }

    std::istream& operator>> (std::istream& s_in, Indeter& X)
    {
        return s_in>>X.name;
    }

} // Givaro
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
