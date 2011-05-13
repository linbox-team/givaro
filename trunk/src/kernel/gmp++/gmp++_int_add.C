// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int_add.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: M. Samama, T. Gautier
// $Id: gmp++_int_add.C,v 1.5 2010-12-20 12:09:37 bboyer Exp $
// ==========================================================================

#include "gmp++/gmp++.h"

namespace Givaro {
//-------------------------------------------------- operator +
Integer& Integer::addin(Integer& res, const Integer& n)
{
  if (isZero(n)) return res;
  if (isZero(res)) return res = n;
  mpz_add( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&res.gmp_rep, (mpz_srcptr)&n.gmp_rep );
  return res;
}
Integer& Integer::addin(Integer& res, const long n)
{
  if (isZero(n)) return res;
  if (isZero(res)) return res = n;
  int sgn = GMP__SGN(n);
  if (sgn >0) mpz_add_ui( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, n);
  else mpz_sub_ui((mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, -n);
  return res;
}
Integer& Integer::addin(Integer& res, const unsigned long n)
{
  if (isZero(n)) return res;
  if (isZero(res)) return res = n;
  mpz_add_ui( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&res.gmp_rep, n);
  return res;
}

Integer& Integer::add(Integer& res, const Integer& n1, const Integer& n2)
{
  if (isZero(n1)) return res = n2;
  if (isZero(n2)) return res = n1;
  mpz_add( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&n1.gmp_rep, (mpz_srcptr)&n2.gmp_rep);
  return res;
}
Integer& Integer::add(Integer& res, const Integer& n1, const long n2)
{
  if (isZero(n1)) return res = n2;
  if (isZero(n2)) return res = n1;
  int sgn = GMP__SGN(n2);
  if (sgn >0) mpz_add_ui( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&n1.gmp_rep, n2);
  else mpz_sub_ui((mpz_ptr)&res.gmp_rep, (mpz_srcptr)&n1.gmp_rep, -n2);
  return res;
}
Integer& Integer::add(Integer& res, const Integer& n1, const unsigned long n2)
{
  if (isZero(n1)) return res = n2;
  if (isZero(n2)) return res = n1;
  mpz_add_ui( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&n1.gmp_rep, n2);
  return res;
}


Integer& Integer::operator += (const Integer& n)
{
  if (isZero(n)) return *this;
  if (isZero(*this)) return logcpy(n);
  mpz_add( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, (mpz_srcptr)&n.gmp_rep) ;
  return *this;
}

Integer& Integer::operator += (const unsigned long l)
{
  if (l==0) return *this;
  if (isZero(*this)) return logcpy(Integer(l));
  mpz_add_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, l);
  return *this;
}

Integer& Integer::operator += (const long l)
{
  if (l==0) return *this;
  if (isZero(*this)) return logcpy(Integer(l));
  int sgn = GMP__SGN(l);
  if (sgn >0) mpz_add_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, l);
  else mpz_sub_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, -l);
  return *this;
}


Integer Integer::operator + (const Integer& n) const
{
  if (isZero(n)) return *this;
  if (isZero(*this)) return n;
  Integer res;
  mpz_add( (mpz_ptr)&(res.gmp_rep), (mpz_srcptr)&gmp_rep, (mpz_srcptr)&n.gmp_rep) ;
  return res;
}

Integer Integer::operator + (const unsigned long l) const
{
  if (l==0) return *this;
  if (isZero(*this)) return Integer(l);
  Integer res;
  mpz_add_ui( (mpz_ptr)&(res.gmp_rep), (mpz_srcptr)&gmp_rep, l);
  return res;
}

Integer Integer::operator + (const long l) const
{
  if (l==0) return *this;
  if (isZero(*this)) return Integer(l);
  Integer res;
  int sgn = GMP__SGN(l);
  if (sgn >0) mpz_add_ui( (mpz_ptr)&(res.gmp_rep), (mpz_srcptr)&gmp_rep, l);
  else mpz_sub_ui( (mpz_ptr)&(res.gmp_rep), (mpz_srcptr)&gmp_rep, -l);
  return res;
}

}
