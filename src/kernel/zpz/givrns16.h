#ifndef _ARITHMODU16_H
#define _ARITHMODU16_H
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givrns16.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givrns16.h,v 1.2 2004-10-11 12:29:50 jgdumas Exp $
// ==========================================================================
// Description:
//  Modular arithmetic for GIVARO. Here is defined arithmetic functions
//  on rns representation and interface between RNS (< 2^16) and Integer.

#include "givaro/givconfig.h"
#include "givaro/giverror.h"
#include "givaro/givarray0.h"
#include "givaro/givzpz16std.h"
#include "givaro/givinteger.h"

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


class RNS16system  {
public:
  typedef Integer   			    Ring;
  typedef ZpzDom<Std16>::Element	    Modulo;
  typedef Array0<Modulo> 		    array;

  // Default Cstor, Dstor/Cstor of recopy: 
  RNS16system();
  ~RNS16system(); 
  RNS16system(const RNS16system& R); 

  // -- Cstor with given primes 
  RNS16system( const array& primes );

  // -- Computation of a mixed-radix representation of the residus.
  void RnsToMixedRadix(array&  mixrad, const array&  residu) const; 

  // -- Convert a mixed radix representation to an Integer
  void MixedRadixToRing( Integer& res,  const array& mixrad ) const;

  // -- Convert an Ring element to a its RNS16 representation
  // with the "this" rns system.
  void RingToRns( array& residu, const Integer& a ) const;

  // -- Fast conversion: requires pre-computation (first time it was called)
  void fastRingToRns( array& residu, const Integer& a ) const;

  // -- Convert a RNS16 representation to an Integer element
  void RnsToRing( Integer& a, const array& residu ) const;

  // -- Fast conversion: requires pre-computation (first time it was called)
  void fastRnsToRing( Integer& a, const array& residu ) const;

  // ------------- Access methods
 
  // -- Returns the number of primes of this ctxt
  int NumOfPrimes() const { return _primes.size(); } 

  // -- Returns a array to the begin of the array of primes
  const array& Primes() const;
  // -- Returns the ith primes of the rns system
  const Modulo ith(const size_t i) const;

  // -- Returns a array of the reciprocal ck = (\prod_{j=0..k-1)p_j)^(-1) [pk]
  const array& Reciprocals() const;
  const Modulo reciprocal(const size_t i) const;

protected:
  // -- Compute some fields of the structure :
  void ComputeCk();

  // -- Compute the Qk for Ring -> RNS, allocate U
  void ComputeQk();

  array  _primes; 	// - array of the relatively primes numbers
  array  _ck;     	// - reciprocals, _ck[0] = 1, same size as _primes 

  // -- for fast conversion
  size_t _sizek;
  size_t _log2k;
  Array0<Integer> _u;	// - cf algo Aho, Ullman & Hopcroft
  Array0<Integer> _qk;	// - cf algo Aho, Ullman & Hopcroft
};

#endif
