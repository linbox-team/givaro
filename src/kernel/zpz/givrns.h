#ifndef _ARITHMODU_H
#define _ARITHMODU_H
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givrns.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givrns.h,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description:
//  Modular arithmetic for GIVARO. Here is defined arithmetic functions
//  on rns representation and interface between RNS and Integer.

#include "givaro/givconfig.h"
#include "givaro/giverror.h"

  // ---------------------------------------------  class RNSsystem
  // Structure which manages list of primes in order to do 
  // convertion between integer and RNS number system using
  // a mixed radix form.
  // This class is parametric by the type of Ring and of Modulo 
  // Modulo type should have following members:
  // - Modulo::zero: the zero of Modulo
  // - Modulo::mul(p, a,b): returns the product a*b [p]  
  // - Modulo::muladd(p, a,b,c): returns the product a*b +c [p] 
  // - Modulo::inv(p, a): returns the multiplicative inverse a^(-1) [p]  
  // - operator= : assignment
  // Ring should have:
  // - Ring( Modulo ): conversion
  // - operator *= (Ring&, const Modulo&)
  // - operator += (Ring&, const Modulo&)
  // - Ring::mod(a, p): returns the (residu) a [p]


template<class RING, class MODULO>
class RNSsystem  {
public:
  typedef RING   ring;
  typedef MODULO modulo;
  typedef Array0<MODULO> array;

  // Default Cstor, Dstor/Cstor of recopy: 
  RNSsystem();
  ~RNSsystem(); 
  RNSsystem(const RNSsystem& R); 

  // -- Cstor with given primes 
  RNSsystem( const array& primes );

  // -- Computation of a mixed-radix representation of the residus.
  void RnsToMixedRadix(array&  mixrad, const array&  residu) const; 

  // -- Convert a mixed radix representation to an Integer
  void MixedRadixToRing( RING& res,  const array& mixrad ) const;

  // -- Convert an Ring element to a its RNS representation
  // with the "this" rns system.
  void RingToRns( array& rns, const RING& a ) const;

  // -- Convert a RNS representation to an RING element
  void RnsToRing( RING& a, const array& rns ) const;

  // ------------- Access methods
 
  // -- Returns the number of primes of this ctxt
  int NumOfPrimes() const { return _primes.size(); } 

  // -- Returns a array to the begin of the array of primes
  const array& Primes() const;
  // -- Returns the ith primes of the rns system
  const MODULO ith(const size_t i) const;

  // -- Returns a array of the reciprocal ck = (\prod_{j=0..k-1)p_j)^(-1) [pk]
  const array& Reciprocals() const;
  const MODULO reciprocal(const size_t i) const;

protected:
  // -------------- Compute some fields of the structure :
  void ComputeCk();

  array  _primes; 	// - array of the primes
  array  _ck;     	// - reciprocals, _ck[0] = 1, same size as _primes 
};

//#include "givaro/givrnscstor.inl"
//#include "givaro/givrnsconvert.inl"

#endif
