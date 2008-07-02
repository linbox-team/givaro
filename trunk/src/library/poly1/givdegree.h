#ifndef _POLY1DEBGREE_H_
#define _POLY1DEBGREE_H_
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givdegree.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givdegree.h,v 1.3 2008-07-02 12:53:02 jgdumas Exp $
// Description: opaque class for Degree of polynomial. Degree of polynomial
// 0 is Degree::deginfty with value DEGPOLYZERO.
// ==========================================================================

#include <iostream>

// -- Degree type for polynomials
// 
class Degree { 
public:
    typedef long value_type;

  enum { DEGPOLYZERO =-1};
    Degree(long a = DEGPOLYZERO): _deg(a) { }
    
// JGD 03.06.2003
// Commented out because of ambiguous overload on the operators
//   // -- cast --> long
//  operator long() { return _deg; }            
//  operator int() { return _deg; }            

  // -- Degree of zero polynomial
  static const long deginfty = DEGPOLYZERO;

  // -- cvrt
  long value() const { return _deg; }

  // -- Basic arithmetic:
  Degree operator+( const Degree& d) const { return Degree(_deg+d._deg); }
  Degree operator-( const Degree& d) const { return Degree(_deg-d._deg); }
  Degree operator*( const Degree& d) const { return Degree(_deg*d._deg); }
  Degree operator/( const Degree& d) const { return Degree(_deg/d._deg); }
  long operator++() { return ++_deg; }
  long operator--() { return --_deg; }
  long operator++(int) { return _deg++; }
  long operator--(int) { return _deg--; }
  
  // -- Comparizon:
  int operator==( const Degree& d) const { return _deg == d._deg; }
  int operator!=( const Degree& d) const { return _deg != d._deg; }
  int operator<=( const Degree& d) const { return _deg <= d._deg; }
  int operator< ( const Degree& d) const { return _deg <  d._deg; }
  int operator>=( const Degree& d) const { return _deg >= d._deg; }
  int operator> ( const Degree& d) const { return _deg >  d._deg; }
  
  // -- methods
    friend std::ostream& operator<< (std::ostream& o, const Degree& d) { return o << d._deg; }
    friend std::istream& operator>> (std::istream& i, Degree& d) { return i >> d._deg; }
    

public:
  long _deg;
};

inline long value(const Degree& d) { return d._deg; }

#endif
