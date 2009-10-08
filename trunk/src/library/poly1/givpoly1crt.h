// ==========================================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software. 
// see the COPYRIGHT file for more details.
// Authors: J-G Dumas
// Time-stamp: <08 Oct 09 12:39:08 Jean-Guillaume.Dumas@imag.fr> 
// Description: Polynomial Chinese Remaindering of degree 1
// ==========================================================================
#ifndef _GIVARO_Poly1_CRT_H
#define _GIVARO_Poly1_CRT_H
#include <givaro/givpoly1.h>
#include <givaro/givindeter.h>

template<class Field>
class Poly1CRT  {
    typedef Poly1CRT<Field> 			Self_t;
public:
    typedef Field				Field_t;
    typedef Poly1Dom<Field, Dense>		Ring_t;
    typedef typename Field_t::Element 		Type_t;
    typedef typename Ring_t::Element 		Element;
    typedef Array0<Type_t> 			array_T;
    typedef Array0<Element> 			array_E;

        // Default Cstor, Dstor/Cstor of recopy: 
    Poly1CRT() ;
    ~Poly1CRT(); 
    Poly1CRT( const Self_t& R); 
    
        // -- Cstor with given residues so that irreds are (X-primes[i])
    Poly1CRT( const Field& F, const array_T& primes, const Indeter& X = Indeter() );
    
        // -- Convert a Ring Element to a its RNS representation
        // with the "this" rns system.
    array_T& RingToRns( array_T& rns, const Element& a ) const;
    
        // -- Convert a RNS representation to a RING Element
    Element& RnsToRing( Element& a, const array_T& rns );

        // ------------- Access methods
 
        // -- Returns the number of primes of this ctxt
    int size() const { return _primes.size(); } 

        // -- Returns a array to the begin of the array of primes
    const array_T& Primes() const;
        // -- Returns the ith primes of the rns system
    const Type_t& ith(const size_t i) const;

        // -- Returns an array of the reciprocal ck = (\prod_{j=0..k-1)p_j)^(-1) [pk]
    const array_E& Reciprocals() const;
    const Element& reciprocal(const size_t i) const;

    const Field_t& getdomain() { return _F; }
    const Ring_t&  getpolydom() { return _PolRing; }
protected:
        // -------------- Compute some fields of the structure :
    void ComputeCk();

    const Indeter  _XIndet;
    const Field_t& _F;
    const Ring_t   _PolRing;
    array_T   _primes; 	// - array of the primes
    array_E       _ck;  // - Radix list reciprocals
};

#include "givaro/givpoly1crtcstor.inl"
#include "givaro/givpoly1crtconvert.inl"

#endif
