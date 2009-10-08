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
#ifndef __GIVARO__PRIMALITY_INL
#define __GIVARO__PRIMALITY_INL
#include <math.h>
#include "givaro/givintprime.h"

// =================================================================== //
// Primality tests and factorization algorithms
// =================================================================== //

// =================================================================== //
// Primality tests 
// =================================================================== //
template<class RandIter> unsigned int IntPrimeDom::Miller(RandIter& g, const Integer& n) const
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

    
template<class RandIter>
IntPrimeDom::Rep& IntPrimeDom::test_Lehmann(RandIter& g, Rep& r, const Rep& n) const {
        // Monte Carlo algorithm
        // returns n-1  : n prime with probability 1/2
        // returns 1    : n composite with probability 1/2
        // else         : n composite
    IntPrimeDom::Rep A;
    random(g,A,n);
    return powmod(r,A,(n-1)/2,n);
}

template<class RandIter>
int IntPrimeDom::Lehmann(RandIter& g, const Rep& n)  const 
{
    if (n < 2) return 0;
    if (n <= 3) return 1;
    IntPrimeDom::Rep tmp;
    IntPrimeDom::test_Lehmann(g,tmp,n);
    if (tmp == (n-1))
        return 1;
    return 0;
}
#endif