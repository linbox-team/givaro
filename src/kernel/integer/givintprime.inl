// =================================================================== //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Givaro : Prime numbers
//              Primality tests
// Time-stamp: <29 Jun 05 14:11:07 Jean-Guillaume.Dumas@imag.fr>
// =================================================================== //
#ifndef __GIVARO_primality_INL
#define __GIVARO_primality_INL
#include <cmath>
#include "givaro/givintprime.h"

// =================================================================== //
// Primality tests and factorization algorithms
// =================================================================== //

// =================================================================== //
// Primality tests
// =================================================================== //

namespace Givaro {

    template<class MyRandIter> unsigned int IntPrimeDom::Miller(MyRandIter& g, const Integer& n) const
    {
            // Monte Carlo algorithm
            // returns 1    : n prime with probability 3/4
            // returns 0    : n composite
        if (n < 2) return 0;
        if (n <= 3) return 1;
        IntPrimeDom::Rep t=n-1,a,q;
        random(g,a,n);
        long s=0;
        for( ; !( (int)t & 0x1) ; t>>=1, ++s) { }
        powmod(q,a,t,n);
        if ( (q==1) || (q == (n-1))) return 1;
            // for(;s>1;--s) {
        for(;--s>0;) {
            q = (q*q) % n;
            if (q == (n-1)) return 1;
        }
        return 0;
    }


    template<class MyRandIter>
    IntPrimeDom::Rep& IntPrimeDom::test_Lehmann(MyRandIter& g, Rep& r, const Rep& n) const {
            // Monte Carlo algorithm
            // returns n-1  : n prime with probability 1/2
            // returns 1    : n composite with probability 1/2
            // else         : n composite
        IntPrimeDom::Rep A;
        random(g,A,n);
        return powmod(r,A,(n-1)/2,n);
    }

    template<class MyRandIter>
    int IntPrimeDom::Lehmann(MyRandIter& g, const Rep& n)  const
    {
        if (n < 2) return 0;
        if (n <= 3) return 1;
        IntPrimeDom::Rep tmp;
        IntPrimeDom::test_Lehmann(g,tmp,n);
        if (tmp == (n-1))
            return 1;
        return 0;
    }


} // Givaro
#endif // __GIVARO_primality_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
