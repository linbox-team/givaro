// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/rational/givratmisc.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: M. Samama
// $Id: givratmisc.C,v 1.5 2009-10-01 09:07:36 jgdumas Exp $
// ==========================================================================
// Description:

#include "givaro/givrational.h"

namespace Givaro {

//
Rational& Rational::reduce()
{
  Integer t = gcd(num, den);
  if (!isOne(t) )
  {
	  num /= t;
	  den /= t;
  }
  return *this;
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

const Integer round(const Rational& x)  // GV Jeu  2 aoû 2018 16:38:58 CEST
{
  Integer q;
  Integer r;
  Integer::divmod(q, r, abs(x.num), abs(x.den));
  r <<= 1; 
  if (absCompare(r,x.den) > 0) // GV was < 0, and changed with abs
      q += 1;
  return sign(x)*q;

}

const Rational pow (const Rational& x, const int64_t y)
{
  Rational r;
  if (y >= 0)
  {
    r.num = pow(x.num, (int64_t) y);
    r.den = pow(x.den, (int64_t) y);
  }
  else
  {
    r.den = pow(x.num, (int64_t) -y);
    r.num = pow(x.den, (int64_t) -y);
    if (sign(r.den) < 0)
    {
      r.num = -r.num ;
      r.den = -r.den ;
    }
  }
  return r;
}

} // namespace Givaro
