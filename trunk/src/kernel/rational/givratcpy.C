// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/rational/givratcpy.C,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: M. Samama
// $Id: givratcpy.C,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description:

#include "givaro/givrational.h"

Rational& Rational::logcpy (const Rational &r)
{ 
  if (this == &r) return *this ;
  num.logcpy(r.num) ; den.logcpy(r.den) ; 
  return *this ;
} 

// same that Rational::logcpy function
Rational& Rational::operator= (const Rational &r)
{
  if (this == &r) return *this ;
  num.logcpy(r.num) ; den.logcpy(r.den) ; 
  return *this ;
} 

Rational& Rational::copy (const Rational &r)
{
  if (this == &r) return *this ;
  num.copy(r.num) ; den.copy(r.den) ; 
  return *this ;
}


 
