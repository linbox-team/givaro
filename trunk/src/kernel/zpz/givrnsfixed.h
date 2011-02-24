// Copyright(c)'1994-2011 by The Givaro group 
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <24 Feb 11 17:20:59 Jean-Guillaume.Dumas@imag.fr>
// ==========================================================================
// Description:
//  Chinese Remainder Algorithm.
#ifndef __GIVARO_arithmodu_fixedprimes_H
#define __GIVARO_arithmodu_fixedprimes_H

#include <vector>
#include "givaro/givrns.h"
#include "givaro/givrandom.h"
#include "givaro/givintprime.h"
#include "givaro/givzpzInt.h"


template<class Ints>
class RNSsystemFixed  {
    typedef RNSsystemFixed<Ints> Self_t;
    typedef RNSsystem<Ints, ZpzDom<Ints> > RNS_t;
public:
    typedef std::vector<Ints>     array;
    typedef std::vector<array>    tree;

        // Default Cstor, Dstor/Cstor of recopy:
    RNSsystemFixed() ;
    ~RNSsystemFixed();
    RNSsystemFixed(const Self_t& R);

        // -- Cstor with given primes
    RNSsystemFixed( const array& primes );

        // -- Convert a RNS representation to a Ints Element
    Ints& RnsToRing( Ints& a, const array& rns ) const;

        // ------------- Access methods

        // -- Returns the number of primes of this ctxt
    int size() const { return _primes.size(); }

        // -- Returns a array to the begin of the array of primes
    const tree& Primes() const;
        // -- Returns the ith primes of the rns system
    const Ints ith(const size_t i) const;


protected:
    Ints& RnsToRingLeft( Ints& I, const array& residues, const int level, const int col ) const;
    Ints& RnsToRingRight( Ints& I, const array& residues, const int level, const int col ) const;

    tree  _primes; 	// - array of the primes and reciprocals
    RNS_t _RNS;		// - unbalanced recovery
};

#include "givaro/givrnsfixed.inl"

#endif
