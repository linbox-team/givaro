// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givrns.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givrns.h,v 1.5 2011-02-02 16:23:56 bboyer Exp $
// ==========================================================================
// Description:
//  Modular arithmetic for GIVARO. Here is defined arithmetic functions
//  on rns representation and interface between RNS and Integer,
//  all is done via the Chinese Remainder Algorithm.
#ifndef __GIVARO_arithmodu_H
#define __GIVARO_arithmodu_H

#include "givaro/givconfig.h"
#include "givaro/giverror.h"
#include "givaro/givarray0.h"

    // ---------------------------------------------  class RNSsystem
    // Structure which manages list of domains in order to
    // convert integer to/from RNS number system using
    // a mixed radix form.
    // This class is parameterized by the type of Ring and of Domain
    // Ring should have:
    // - Ring( Modulo ) and Domain.init(Ring) conversions
    // - operator *= (Ring&, const Modulo&)
    // - operator += (Ring&, const Modulo&)

template<class RING, class Domain>
class RNSsystem  {
    typedef RNSsystem<RING, Domain> Self_t;
public:
    typedef RING   ring;
    typedef typename Domain::Element modulo;
    typedef Array0<modulo> array;
    typedef Array0<Domain> domains;

        // Default Cstor, Dstor/Cstor of recopy:
    RNSsystem() ;
    ~RNSsystem() ;
    RNSsystem(const Self_t& R);

        // -- Cstor with given primes
    RNSsystem( const domains& primes );

        // -- Computation of a mixed-radix representation of the residus.
    void RnsToMixedRadix(array&  mixrad, const array&  residu) const;

        // -- Convert a mixed radix representation to an Integer
    RING& MixedRadixToRing( RING& res,  const array& mixrad ) const;

        // -- Convert a Ring Element to a its RNS representation
        // with the "this" rns system.
    void RingToRns( array& rns, const RING& a ) const;

        // -- Convert a RNS representation to a RING Element
    RING& RnsToRing( RING& a, const array& rns ) const;

        // ------------- Access methods

        // -- Returns the number of primes of this ctxt
    int size() const { return _primes.size(); }

        // -- Returns a array to the begin of the array of primes
    const domains& Primes() const;
        // -- Returns the ith primes of the rns system
    const Domain ith(const size_t i) const;

        // -- Returns an array of the reciprocal ck = (\prod_{j=0..k-1)p_j)^(-1) [pk]
    const array& Reciprocals() const;
    const modulo reciprocal(const size_t i) const;

        // -- Cstor with given primes
    void setPrimes( const domains& primes );
protected:
        // -------------- Compute some fields of the structure :
    void ComputeCk();

    domains  _primes; 	// - array of the primes
    array  _ck;     	// - reciprocals, _ck[0] = 1, same size as _primes
};

#include "givaro/givrnscstor.inl"
#include "givaro/givrnsconvert.inl"

#endif // __GIVARO_arithmodu_H