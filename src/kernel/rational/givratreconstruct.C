// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/rational/givratreconstruct.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: Jean-Guillaume Dumas
// $Id: givratreconstruct.C,v 1.4 2010-04-12 15:54:39 jgdumas Exp $
// ==========================================================================
// Description:

#include "givaro/givrational.h"
#include <iostream>

namespace Givaro {


    Rational::Rational(const Integer& f, const Integer& m, const Integer& k, bool recurs )
    {
        bool res = this->ratrecon(f,m,k,recurs);
        if (recurs)
            for( Integer newk = k + 1; (!res) && (newk<f) ; ++newk)
                res = this->ratrecon(f,m,newk,true);
    }



    bool Rational::ratrecon(const Integer& f, const Integer& m, const Integer& k, bool recurs )
    {

        return Givaro::ratrecon(this->num,this->den,f,m,k,true,recurs);
    }


} // namespace Givaro

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
