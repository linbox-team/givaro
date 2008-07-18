// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/rational/givratmisc.C,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: M. Samama
// $Id: givratmisc.C,v 1.3 2008-07-18 12:42:37 jgdumas Exp $
// ==========================================================================
// Description:

#include "givaro/givrational.h"

// 
void Rational::reduce(void) 
{
  Integer t = gcd(num, den);
  if (!isOne(t) )
  {
	  num /= t;
	  den /= t;
   }
}
	  
const Integer trunc(const Rational &r)
{
  return r.num / r.den;
}

const Integer floor(const Rational &x)
{
  Integer q, r;
  Integer::divmod(q, r, x.num, x.den);
  if (sign(x.num) < 0 && sign(r) != 0) q -= 1;
  return q;
}

const Integer ceil(const Rational &x)
{
  Integer r, q;

  Integer::divmod(q, r, x.num, x.den);
  if (sign(x.num) >= 0 && sign(r) != 0) q += 1;
  return q;
}

const Integer round(const Rational& x) 
{
  Integer q;
  Integer r;
  Integer::divmod(q, r, x.num, x.den);
  r <<= 1UL;
  if (absCompare(r,x.den) < 1)
  {
    if (sign(x.num) >= 0)
      q += 1;
    else
      q -= 1;
  }
  return q;
}

const Rational pow (const Rational& x, const long y)
{
  Rational r;
  if (y >= 0)
  {
    r.num = pow(x.num, (long) y);
    r.den = pow(x.den, (long) y);
  }
  else
  {
    r.den = pow(x.num, (long) -y);
    r.num = pow(x.den, (long) -y);
    if (sign(r.den) < 0)
    {
      r.num = -r.num ; 
      r.den = -r.den ;
    }
  }
  return r;
}

