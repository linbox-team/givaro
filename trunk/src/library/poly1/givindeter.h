#ifndef _INDETER_H_
#define _INDETER_H_
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givindeter.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givindeter.h,v 1.2 2008-02-21 18:42:48 pernet Exp $
// ==========================================================================
// Description:
// - indeterminates for polynomial manipulation

#include <iostream>

#include "givaro/givarray0.h"

class Indeter {
public :
     
  // -- Cstor: recopy the string
  Indeter(const char *s = "");
  // -- Cstor of recopy
  Indeter(const Indeter& s) : name(s.name, givNoCopy()){};

  // -- Dstor
  virtual ~Indeter();
 
  // -- assignement
  Indeter& operator= (const Indeter& s);

  // -- Comparizon operators:
  // all comparizons are based on this virtual method,
  // which returns : -1 iff *this < b, 0 iff *this == b and
  // +1 else. This comparizon method gives the natural order
  // for multivariate polynomials.
virtual int compare(const Indeter& b)  const;

  // -- methods	
  friend std::ostream& operator<< (std::ostream& o, const Indeter& X);
  friend std::istream& operator>> (std::istream& o, Indeter& X);

protected:
  Array0<char> name;
};

  // Inline members functions :
inline int operator==(const Indeter& i1, const Indeter &i2) 
  { return i1.compare(i2) ==0; }

inline int operator!=(const Indeter& i1, const Indeter &i2) 
  { return i1.compare(i2) !=0; }

inline int operator<= (const Indeter& i1, const Indeter &i2)  
  { return i1.compare(i2) <=0; }

inline int operator<  (const Indeter& i1, const Indeter &i2) 
  { return i1.compare(i2) <0; }

inline int operator>= (const Indeter& i1, const Indeter &i2) 
  { return i1.compare(i2) >=0; }

inline int operator>  (const Indeter& i1, const Indeter &i2)
  { return i1.compare(i2) >0; }

#endif
