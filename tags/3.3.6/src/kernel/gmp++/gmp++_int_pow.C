// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int_pow.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: JG Dumas
// $Id: gmp++_int_pow.C,v 1.4 2009-09-17 14:28:22 jgdumas Exp $
// ==========================================================================
// Description:

#include "gmp++/gmp++.h"
int isperfectpower(const Integer& n)
{
	return mpz_perfect_power_p((mpz_ptr)&(n.gmp_rep));
}

Integer& pow(Integer& Res, const unsigned long n, const unsigned long p)
{
  mpz_ui_pow_ui( (mpz_ptr)&(Res.gmp_rep), n, p);
  return Res;
}
Integer& pow(Integer& Res, const Integer& n, const unsigned long p)
{
__gmpz_pow_ui( (mpz_ptr)&(Res.gmp_rep), (mpz_ptr)&n.gmp_rep, p);
  return Res;
}

Integer pow(const Integer& n, const unsigned long p)
{
  if (p == 0) return Integer::one;

  Integer Res;
  return pow(Res,n,p);
}

Integer& pow(Integer& Res, const Integer& n, const long l)
{
        // Beware of negative values
	return pow(Res, n, (unsigned long) GMP__ABS(l) );
}
Integer pow(const Integer& n, const long l)
{
  if (l < 0)  return Integer::zero;
  return pow(n, (unsigned long) GMP__ABS(l) );
}

Integer& powmod(Integer& Res, const Integer& n, const unsigned long p, const Integer& m)
{
  mpz_powm_ui( (mpz_ptr)&(Res.gmp_rep), (mpz_ptr)&n.gmp_rep, p, (mpz_ptr)&m.gmp_rep);
  return Res;
}

Integer powmod(const Integer& n, const unsigned long p, const Integer& m)
{
  if (p == 0) return Integer::one;
  Integer Res;
  return powmod(Res,n,p,m);
}

Integer& powmod(Integer& Res, const Integer& n, const long e, const Integer& m)
{
    if (e < 0) {
        inv(Res, n, m);
        return powmod(Res, Res, (unsigned long)GMP__ABS(e), m);
    } else
        return powmod (Res, n, (unsigned long)(e), m);
}
Integer powmod(const Integer& n, const long e, const Integer& m)
{
    Integer Res;
    return powmod(Res, n, e, m);
}


Integer& powmod(Integer& Res, const Integer& n, const Integer& e, const Integer& m)
{
  mpz_powm( (mpz_ptr)&(Res.gmp_rep), (mpz_ptr)&n.gmp_rep, (mpz_ptr)&e.gmp_rep, (mpz_ptr)&m.gmp_rep);
  return Res;
}
Integer powmod(const Integer& n, const Integer& e, const Integer& m)
{
  if (e == 0) return Integer::one;
  if (e < 0)  return Integer::zero;
  Integer Res;
  return powmod(Res, n, e, m);
}
