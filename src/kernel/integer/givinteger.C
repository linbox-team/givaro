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

namespace Givaro {

    void ZRing<Integer>::RationalReconstruction(Rep& a, Rep& b, const Rep& f, const Rep& m, const Rep& k, bool forcereduce, bool recurs) const {
        Rational::RationalReconstruction(a,b,f,m,k,forcereduce,recurs);
    }
    bool ZRing<Integer>::ratrecon(Rep& num, Rep& den, const Rep& f, const Rep& m, const Rep& k, bool forcereduce, bool recurs) const {
        return Rational::ratrecon(num,den,f,m,k,forcereduce,recurs);
    }
    void ZRing<Integer>::reconstructRational (Element& a, Element& b, const Element& x, const Element& m) const
    {this->RationalReconstruction(a,b, x, m, Givaro::sqrt(m), true, true);}
    void ZRing<Integer>::reconstructRational (Element& a, Element& b, const Element& x, const Element& m, const Element& bound) const
    {this->RationalReconstruction(a,b, x, m, bound, true, true);}
    bool ZRing<Integer>::reconstructRational (Element& a, Element& b, const Element& x, const Element& m, const Element& a_bound, const Element& b_bound) const
    {
        Element bound = x/b_bound;
        this->RationalReconstruction(a,b,x,m, (bound>a_bound?bound:a_bound), true, false);
        return b <= b_bound;
    }

}
