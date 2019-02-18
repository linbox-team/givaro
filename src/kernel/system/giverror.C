// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/system/giverror.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: giverror.C,v 1.3 2009-09-17 14:28:23 jgdumas Exp $
// ==========================================================================
// Description:
// - error exception

#include "givaro/giverror.h"
#include <iostream>

namespace Givaro {

    std::ostream& GivError::print( std::ostream& o ) const
    { return o << strg ; }


    GivError::~GivError(){}

    GivMathError::~GivMathError(){}

    GivBadFormat::~GivBadFormat(){}

    GivMathDivZero::~GivMathDivZero(){}

    void GivError::throw_error( const GivError& err )
    {
        throw err;
    }

    std::ostream& operator<< (std::ostream& o, const GivError& E)
    {
        return E.print(o) ;
    }
} // namespace Givaro
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
