// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/rational/givratmuldiv.C,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: M. Samama
// $Id: givratmuldiv.C,v 1.2 2005-06-14 14:53:14 pernet Exp $
// ==========================================================================
// Description:

#include "givaro/givrational.h"
#include "givaro/giverror.h"

// --------------------------------------- operator * 
Rational Rational::operator * (const Rational& r) const
{

  if (isZero(r)) return Rational(0L) ;
  if (isZero(*this)) return Rational(0L) ;
  if (isOne(r)) return *this ;
  if (isOne(*this)) return r ;
  if (isinteger(*this) && isinteger(r))
    return Rational(num*r.num) ;

  if (absCompare(den, r.den) == 0)
    return Rational(num * r.num, den * r.den, 0) ;

  if (Rational::flags == Rational::NoReduce)
  {
     return Rational( num*r.num, den*r.den, 0) ;
  }

  Integer d1 = gcd(num, r.den);
  Integer d2 = gcd(den, r.num);
  return Rational( (num / d1) * (r.num / d2), (den / d2) * (r.den / d1), 0 );
}


// --------------------------------------- operator *= 
Rational& Rational::operator *= (const Rational& r)
{

  if (isZero(r)) return *this=Rational(0L) ;
  if (isZero(*this)) return *this ;
  if (isOne(r)) return *this ;
  if (isOne(*this)) return *this=r ;
  if (isinteger(*this) && isinteger(r)) {
      num *= r.num;
      return *this;
  }

  if ( (absCompare(den, r.den) == 0) || (Rational::flags == Rational::NoReduce) ) {
      num *= r.num;
      den *= r.den;
      return *this;
  }
  
  Integer d1 = gcd(num, r.den);
  Integer d2 = gcd(den, r.num);

  num /= d1;
  num *= (r.num / d2);
  den /= d2;
  den *= (r.den / d1);
  
  return *this;
}


// --------------------------------------- operator / 
Rational Rational::operator / (const Rational& r) const 
{
  if ( isZero(r) ) {
     throw GivMathDivZero("*** division by zero, in operator / (const Rational&)") ;
  }
  if (isZero(*this)) return Rational(0L) ;
  if (isOne(r)) return *this ;
  if (isOne(*this))  {
    if (sign(r) < 0)
      return Rational(r.den, r.num, 0) ;
    else 
      return Rational(-r.den, -r.num, 0) ;
  }

  if (absCompare(den, r.den) == 0) 
    return Rational(num, r.num) ;

  if (Rational::flags == Rational::NoReduce)
     return Rational( num*r.den, den*r.num, 0) ;
  
  Integer d1 = gcd(num, r.num);
  Integer d2 = gcd(den, r.den);
  Integer resnum = (num / d1) * (r.den / d2);
  if (sign(r.num) < 0)
    resnum = -resnum ;
  Integer resden = (den / d2) * (r.num / d1);
    //rden can't be nul
  if (sign(resden) <0) resden = abs(resden) ;
  return Rational(resnum, resden, 0);
}

  
// --------------------------------------- operator /= 
Rational& Rational::operator /= (const Rational& r)
{
  if ( isZero(r) ) {
     throw GivMathDivZero("*** division by zero, in operator / (const Rational&)") ;
  }
  if (isZero(*this)) return *this ;
  if (isOne(r)) return *this ;
  if (isOne(*this))  {
      if (sign(r) < 0) {
          num = r.den;
          den = r.num;
          return *this;
      } else {
          num = -r.den;
          den = -r.num;
          return *this;
      }
  }
  
  if (absCompare(den, r.den) == 0) {
      den = r.num;
      return *this;
  }

  if (Rational::flags == Rational::NoReduce) {
      num *= r.den;
      den *= r.num;
      return *this;
  }
  
  Integer d1 = gcd(num, r.num);
  Integer d2 = gcd(den, r.den);

  num /= d1;
  num *= (r.den / d2);
  if (sign(r.num) < 0)
    num = -num ;

  den /= d2;
  den *= (r.num / d1);
  
    //rden can't be nul
  if (sign(den) <0) den = abs(den) ;
  return *this;
}

  
