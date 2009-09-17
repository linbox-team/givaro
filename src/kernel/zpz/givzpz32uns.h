// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givzpz32uns.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software. 
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givzpz32uns.h,v 1.8 2009-09-17 14:28:23 jgdumas Exp $
// ==========================================================================
//
//  Modified by Pascal Giorgi on 2002/02/13  (pascal.giorgi@ens-lyon.fr)
//
// Description:
//   Arithmetic on Z/pZ, with p a prime number less than 2^32
//   Modulo typedef is a signed long number. In case it was modified
//   then bezout algorithm must be changed (coefficient can be negative).
#ifndef _GIVARO_ZPZ32UNSIGNED_H_ 
#define _GIVARO_ZPZ32UNSIGNED_H_ 

#include "givaro/givbasictype.h"
#include "givaro/giverror.h"
#include "givaro/givzpz.h"
#include "givaro/giv_randiter.h"
#include <math.h>


// ==========================================================================
// -- This class implement the standard arithmetic with Modulo Elements:
// - The representation of an integer a in Zpz is the value a % p
// ==========================================================================

template<>
class ZpzDom<Unsigned32> {
public:
  // ----- Exported Types and constantes
  typedef Unsigned32::type Residu_t;         // - type to store residue
  enum { size_rep = sizeof(Residu_t) };      // - size of the storage type
  // ----- Representation of Element of the domain ZpzDom
  typedef Unsigned32::type Rep;
  typedef Unsigned32::type Element;

  // ----- Representation of vector of the Element
  typedef Rep* Array;
  typedef const Rep* constArray;

  // ----- Constantes 
  const Rep zero;
  const Rep one;

  // ----- Constructor 
  ZpzDom();
  ZpzDom( Residu_t p );
  ZpzDom( const ZpzDom<Unsigned32>& F);

  int operator==( const ZpzDom<Unsigned32>& BC) const { return _p == BC._p;}
  int operator!=( const ZpzDom<Unsigned32>& BC) const { return _p != BC._p;}
 
  ZpzDom<Unsigned32>& operator=( const ZpzDom<Unsigned32>& F) { 
      this->_p = F._p; 
      this->_dp = F._dp;
      return *this;
  }

  // ----- Access to the modulus 
  Residu_t residu() const;
  Residu_t size() const {return _p;}
  Rep access( const Rep a ) const { return a; }
  Residu_t characteristic() const { return _p; }
  Integer& characteristic(Integer& p) const { return p=_p; }
  Residu_t cardinality() const { return _p; }


  // ----- Access to the modulus 
  Rep& init( Rep& a ) const;
  void init( const size_t, Array a, constArray b ) const;
  Rep& init( Rep& r , const long a) const ;
  Rep& init( Rep& r , const unsigned long a) const ;
  Rep& init( Rep& a, const int i) const ;
  Rep& init( Rep& a, const unsigned int i) const ;
    Rep& init( Rep& a, const Integer& i) const ;


  // Initialisation from double ( added for FFLAS usage) (C Pernet)
  Rep& init( Rep& a, const double i) const;
  Rep& init( Rep& a, const float i) const;
  // Conversion to double ( added for FFLAS usage) (C Pernet)
  float& convert(float& r, const Rep a ) const { return r = (float)a ;}
  double& convert(double& r, const Rep a ) const { return r = (double)a ;}
  long int& convert(long int& r, const Rep a) const { return r = (long int)a;}
  unsigned long int& convert(unsigned long int& r, const Rep a) const { return r = (unsigned long int)a;}
  int32& convert(int32& r, const Rep a) const { return r = (int32)a;}
  uint32& convert(uint32& r, const Rep a) const { return r = (uint32)a;}
    Integer& convert(Integer& i, const Rep a) const {
        unsigned long ur;
        return i = (Integer)convert(ur, a);
    }        

  // ----- Misc methods
  int isZero( const Rep a ) const;
  int isOne ( const Rep a ) const;
  size_t length ( const Rep a ) const;

  // ----- Equality between two Elements
  int areEqual(const  Rep& a, const Rep& b) const { return a==b;}

  // ----- Operations with reduction: r <- a op b mod p, r <- op a mod p
  Rep& mul (Rep& r, const Rep a, const Rep b) const;
  Rep& div (Rep& r, const Rep a, const Rep b) const;
  Rep& add (Rep& r, const Rep a, const Rep b) const;
  Rep& sub (Rep& r, const Rep a, const Rep b) const;
  Rep& neg (Rep& r, const Rep a) const;
  Rep& inv (Rep& r, const Rep a) const;

  Rep& mulin (Rep& r, const Rep a) const;
  Rep& divin (Rep& r, const Rep a) const;
  Rep& addin (Rep& r, const Rep a) const;
  Rep& subin (Rep& r, const Rep a) const;
  Rep& negin (Rep& r) const;
  Rep& invin (Rep& r) const;

  // ----- Operations with reduction: r <- a op b mod p, r <- op a mod p
  void mul (const size_t sz, Array r, constArray a, constArray b) const;
  void mul (const size_t sz, Array r, constArray a, Rep b) const;

  void div (const size_t sz, Array r, constArray a, constArray b) const;
  void div (const size_t sz, Array r, constArray a, Rep b) const;

  void add (const size_t sz, Array r, constArray a, constArray b) const;
  void add (const size_t sz, Array r, constArray a, Rep b) const;

  void sub (const size_t sz, Array r, constArray a, constArray b) const;
  void sub (const size_t sz, Array r, constArray a, Rep b) const;

  void neg (const size_t sz, Array r, constArray a) const;
  void inv (const size_t sz, Array r, constArray a) const;

  // -- axpy: r <- a * x + y mod p
  Rep& axpy  (Rep& r, const Rep a, const Rep b, const Rep c) const;
  void axpy 
   (const size_t sz, Array r, constArray a, constArray x, constArray c) const;
  // -- axpyin: r <- r + a * x mod p
  Rep& axpyin(Rep& r, const Rep a, const Rep b) const;
  void axpyin 
   (const size_t sz, Array r, constArray a, constArray x) const;

  // -- axmy: r <- a * x - y mod p
  Rep& axmy  (Rep& r, const Rep a, const Rep b, const Rep c) const;
  void axmy 
   (const size_t sz, Array r, constArray a, constArray x, constArray c) const;
  // -- axmyin: r <- r - a * x mod p
  Rep& axmyin(Rep& r, const Rep a, const Rep b) const;
  void axmyin 
   (const size_t sz, Array r, constArray a, constArray x) const;

  // -- Misc: r <- a mod p
  void assign ( const size_t sz, Array r, constArray a ) const;
/* JGD 26.10.99
  void assign ( Rep& r, const Rep a) const;
  void assign ( Rep& r, const long a ) const;
  void assign ( Rep& r, const unsigned long a ) const;
  void assign ( Rep& r, const int a ) const;
  void assign ( Rep& r, const unsigned int a ) const;
*/
  Rep& assign ( Rep& r, const Rep a) const;
  Rep& assign ( Rep& r, const long a ) const;
  Rep& assign ( Rep& r, const unsigned long a ) const;
  Rep& assign ( Rep& r, const short a ) const;
  Rep& assign ( Rep& r, const unsigned short a ) const;
   // ----- random generators
//     Rep& NONZEROGIVRANDOM(Rep&) const ;
//     Rep& GIVRANDOM(Rep&) const ;
    template< class RandIter > Rep& random(RandIter&, Rep& r) const ;
    template< class RandIter > Rep& random(RandIter&, Rep& r, long s) const ;
    template< class RandIter > Rep& random(RandIter&, Rep& r, const Rep& b) const ;
    template< class RandIter > Rep& nonzerorandom(RandIter&, Rep& r) const ;
    template< class RandIter > Rep& nonzerorandom(RandIter&, Rep& r, long s) const ;
    template< class RandIter > Rep& nonzerorandom(RandIter&, Rep& r, const Rep& b) const ;
    
    typedef GIV_randIter< ZpzDom<Unsigned32> , Rep > randIter;

  // <- \sum_i a[i], return 1 if a.size() ==0,
  Rep& reduceadd ( Rep& r, const size_t sz, constArray a ) const; 

  // <- \prod_i a[i], return 1 if a.size() ==0,
  Rep& reducemul ( Rep& r, const size_t sz, constArray a ) const; 

  // <- \sum_i a[i] * b[i] 
  Rep& dotprod ( Rep& r, const size_t sz, constArray a, constArray b ) const; 
  Rep& dotprod ( Rep& r, const int bound, const size_t sz, constArray a, constArray b ) const; 

  // ----- a -> r: uint32 to double
  void i2d ( const size_t sz, double* r, constArray a ) const; 

  // ----- a -> r % p: double to uint32 % p
  void d2i ( const size_t sz, Array r, const double* a ) const; 

  // --- IO methods
  std::istream& read ( std::istream& s );
  std::ostream& write( std::ostream& s ) const;
  std::istream& read ( std::istream& s, Rep& a ) const;
  std::ostream& write( std::ostream& s, const Rep a ) const;

protected:
  // -- based for modular inverse, d = a*u + b*v
//   static const int32 gcdext ( int32& u, int32& v, const int32 a, const int32 b );
  int32& gcdext (int32& d, int32& u, int32& v, const int32 a, const int32 b ) const;
  uint32& invext (uint32& u, const uint32 a, const uint32 b ) const;

protected:
  // -- data representation of the domain:
    Residu_t _p;
    double _dp;

    static void Init();
    static void End();
};


#include "givaro/givzpz32uns.inl"

#endif
