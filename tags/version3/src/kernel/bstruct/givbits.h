#ifndef _BITS_H_
#define _BITS_H_
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/bstruct/givbits.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Author: T. Gautier
// $Id: givbits.h,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description: 
// - field of n bits, for any n

#include <iostream>
#include "givaro/givaromm.h"
#include "givaro/givarray0.h"

class Bits {
public:
  typedef unsigned long base ;

  Bits () ;
  Bits (const size_t n) ; // -- n is the number of bits reclaimed
  Bits( const Bits& B ) ; // -- logical recopy
  ~Bits() ;

  // -- All binary operand must have same type:
  // <op>in are for inplace operator (this receives the result)
  
  const Bits operator&( const Bits& A ) const ;
  Bits& andin( const Bits& A, const Bits B) ;
  Bits& operator&= (const Bits& A) ;
  //inline const Bits operator& (const Bits& A) const { return operator&(A) ; }

  // -- Or 
  const Bits operator|( const Bits& A ) const ;
  Bits& orin( const Bits& A, const Bits B) ;
  Bits& operator|=( const Bits& A )  ;
  //inline const Bits operator| (const Bits& A) const { return operator|(A) ; }

  // -- Or 
  const Bits operator^( const Bits& A ) const ;
  Bits& xorin( const Bits& A, const Bits B) ;
  Bits& operator^=( const Bits& A )  ;
  //inline const Bits operator^ (const Bits& A) const { return operator^(A) ; }

  // -- Not
  const Bits operator~( ) const ;
  Bits& notin( const Bits& A ) ;
  //inline const Bits operator~ () const { return operator~() ; }

  // -- Assignement: physical recopy
  Bits& copy( const Bits& src ) ;
  Bits& operator= (const Bits& B);

  // -- Logical recopy
  Bits& logcopy( const Bits& src ) ;

  // -- Returns the number of non zero bits :
  long numone() const ;

  // -- Returns the index of non zero bits 
  void indexofone( Array0<Bits::base>& ) const ;

  // -- Returns the length (in bit) of this :
  size_t length() const ;
  
  // -- set to 0 all bits:
  void clear() ;
  // -- set to 0 the i-th bit
  void clear(const int i) ;
  // -- set to 1 all bits
  void set() ;
  // -- set to 1 the i-th bit
  void set(const int i) ;
  // -- returns the i-th bit
  int  get(const int i) const ;
  int  operator[] (const int i) const ;

  // -- IO/methods
  std::ostream& print( std::ostream& o ) const ;

protected:
  typedef Array0<Bits::base> Rep ;
  Bits( const Rep& r ) ;
  Rep rep ;

private:
  static GivModule Module;
  friend class GivModule;
  static void Init(int*, char***);  
  static void End();
} ;

#include "givaro/givbits.inl"

#endif
