// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/rational/givratcompare.C,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: M. Samama
// $Id: givratcompare.C,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description:
#include "givaro/givrational.h"

// -------------------------------------------- compare
// returns 1 if a > b, 0 if a == b and -1 otherwise.
int compare (const Rational& a, const Rational &b) 
{ 
  if (iszero(a.num) && iszero(b.num)) 
      return 0 ;
  if (iszero(a.num)) 
      return -sign(b.num);
  if (iszero(b.num))
      return sign(a.num);
  if (sign(a.num) != sign(b.num)) 
      return (sign(a.num) == -1 ? -1 : 1 )  ;
  if (sign(a.num) > 0) 
      return absCompare(a,b) ;
  else return -absCompare(a,b) ;
} 

// -------------------------------------------- absCompare
int absCompare (const Rational& a, const Rational& b)
{ 
  int cnum = absCompare(a.num, b.num) ;
  int cden = absCompare(a.den, b.den) ;
 
  if ( (cnum == -1) && (cden == 1) )
    return -1 ;
  if ( (cnum == 1) && (cden == -1) )
    return 1 ;
  if (cnum == 0)
    return -cden ;
  if (cden == 0)
    return cnum ;

  Integer p1 = a.num * b.den;
  Integer p2 = a.den * b.num;
  return absCompare(p1,p2);
} 

