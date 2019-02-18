// Copyright(c)'1994-2011 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <01 Apr 11 15:43:07 Jean-Guillaume.Dumas@imag.fr>
// ==========================================================================

/*! @file givrnsfixed.h
 * @ingroup zpz
 * @brief Chinese Remainder Algorithm.
 */

#ifndef __GIVARO_arithmodu_fixedprimes_H
#define __GIVARO_arithmodu_fixedprimes_H

#include "givaro/givrns.h"
#include "givaro/givrandom.h"
#include "givaro/givintprime.h"
#include "givaro/modular-integer.h"
#include <vector>

namespace Givaro {


	/*! @brief NO DOC
	 */
template<class Ints>
class RNSsystemFixed  {
    typedef RNSsystemFixed<Ints> Self_t;
    typedef RNSsystem<Ints, Modular<Ints> > RNS_t;
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
    template<class smallIntVector>
    Ints& RnsToRing( Ints& a, const smallIntVector& rns ) ;

        // ------------- Access methods

        // -- Returns the number of primes of this ctxt
    int size() const { return _primes.size(); }

        // -- Returns a array to the beginning of the array of primes
    const tree& Primes() const;
        // -- Returns the ith primes of the rns system
    const Ints ith(const size_t i) const;


protected:
    template<class smallIntVector>
    Ints& RnsToRingLeft( Ints& I, const smallIntVector& residues, const int level, const int col ) ;
    template<class smallIntVector>
    Ints& RnsToRingRight( Ints& I, const smallIntVector& residues, const int level, const int col ) ;

    tree  _primes; 	// - array of the primes and reciprocals
    RNS_t _RNS;		// - unbalanced recovery
};

} // namespace Givaro

#include "givaro/givrnsfixed.inl"

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
