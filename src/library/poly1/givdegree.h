// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givdegree.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givdegree.h,v 1.7 2011-02-02 16:23:56 bboyer Exp $
// ==========================================================================
/** @file givdegree.h
 * @ingroup poly1
 * @brief NO DOC
 * opaque class for Degree of polynomial. Degree of polynomial
 * 0 is Degree::deginfty with value DEGPOLYZERO.
 *
 */

#ifndef __GIVARO_poly1degree_H
#define __GIVARO_poly1degree_H

#include <iostream>

namespace Givaro {
//! Degree type for polynomials
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
  static const long deginfty;

  // -- cvrt
  long value() const { return _deg; }

  // -- Basic arithmetic:
  Degree operator+( const Degree& d) const { return Degree(_deg+d._deg); }
  Degree operator-( const Degree& d) const { return Degree(_deg-d._deg); }
  Degree operator*( const Degree& d) const { return Degree(_deg*d._deg); }
  Degree operator/( const Degree& d) const { return Degree(_deg/d._deg); }
  Degree& operator+=( const Degree& d) { _deg+=d._deg; return *this; }
  Degree& operator-=( const Degree& d) { _deg-=d._deg; return *this; }
  Degree& operator*=( const Degree& d) { _deg*=d._deg; return *this; }
  Degree& operator/=( const Degree& d) { _deg/=d._deg; return *this; }
    Degree operator<<( const int i) const { return Degree(_deg<<i); }
    Degree operator>>( const int i) const { return Degree(_deg>>i); }
    Degree& operator <<=( const int i) { _deg<<=i; return *this;}
    Degree& operator >>=( const int i) { _deg>>=i; return *this;}
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
  int operator==( const long& d) const { return _deg == d; }
  int operator!=( const long& d) const { return _deg != d; }
  int operator<=( const long& d) const { return _deg <= d; }
  int operator< ( const long& d) const { return _deg <  d; }
  int operator>=( const long& d) const { return _deg >= d; }
  int operator> ( const long& d) const { return _deg >  d; }

  // -- methods
    friend std::ostream& operator<< (std::ostream& o, const Degree& d) { return o << d._deg; }
    friend std::istream& operator>> (std::istream& i, Degree& d) { return i >> d._deg; }


public:
  long _deg;
};

//! value
inline long value(const Degree& d) { return d._deg; }
} // Givaro

#endif // __GIVARO_poly1degree_H
