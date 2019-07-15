// ==========================================================================
// $Source: Givaro/src/kernel/integer/givinteger.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: Jean-Guillaume Dumas
// $Id: givratreconstruct.C,v 1.4 2010-04-12 15:54:39 jgdumas Exp $
// ==========================================================================
// Description:

#include "givaro/givinteger.h"
#include "givaro/givrational.h"

namespace Givaro {
    bool ZRing<Integer>::ratrecon(Rep& num, Rep& den, const Rep& f, const Rep& m, const Rep& k, bool forcereduce, bool recurs) const {
        return Rational::ratrecon(num,den,f,m,k,forcereduce,recurs);
    }
    bool ZRing<Integer>::RationalReconstruction
    (Element& a, Element& b, const Element& x, const Element& m) const {
        return Rational::RationalReconstruction(a, b, x, m);
    }
    bool ZRing<Integer>::RationalReconstruction
    (Element& a, Element& b, const Element& x, const Element& m, 
     const Element& bound, bool forcereduce, bool recurs) const {
        return Rational::RationalReconstruction(a,b, x, m, bound, forcereduce, recurs);
    }
    bool ZRing<Integer>::RationalReconstruction
    (Element& a, Element& b, const Element& x, const Element& m, 
     const Element& a_bound, const Element& b_bound) const {
        return Rational::RationalReconstruction(a,b, x, m, a_bound,b_bound);
    }
}

