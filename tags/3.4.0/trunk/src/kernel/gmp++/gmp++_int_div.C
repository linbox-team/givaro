// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int_div.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: M. Samama, T. Gautier
// $Id: gmp++_int_div.C,v 1.9 2011-01-06 18:02:37 bboyer Exp $
// ==========================================================================

#include "gmp++/gmp++.h"

namespace Givaro {

//-------------------------------------------------- operator /
Integer& Integer::divin(Integer& res, const Integer& n)
{
	//  if (isZero(n)) {
	//    GivMathDivZero("[Integer::/]: division by zero");
	//  }
	if (isZero(res)) return res;
	mpz_tdiv_q( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n.gmp_rep );
	return res;
}

Integer& Integer::divin(Integer& res, const long n)
{
	//  if (n ==0) {
	//    GivMathDivZero("[Integer::/]: division by zero");
	//  }
	if (isZero(res)) return res;
	int sgn = GMP__SGN(n);
	mpz_tdiv_q_ui( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, GMP__ABS(n));
	if (sgn <0) return res = -res;
	// if (n<0)
	// mpz_fdiv_q_ui( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&res.gmp_rep, n);
	// else
	// mpz_cdiv_q_ui( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&res.gmp_rep, n);


	return res;
}

Integer& Integer::divin(Integer& res, const unsigned long n)
{
	//  if (n ==0) {
	//    GivMathDivZero("[Integer::/]: division by zero");
	//  }
	if (isZero(res)) return res;
	mpz_tdiv_q_ui( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&res.gmp_rep, n);
	return res;
}

Integer& Integer::div(Integer& res, const Integer& n1, const Integer& n2)
{
	if (isZero(n1)) return res = Integer::zero;
	//  if (isZero(n2)) {
	//    GivMathDivZero("[Integer::/]: division by zero");
	//  }
	mpz_tdiv_q( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, (mpz_ptr)&n2.gmp_rep);
	return res;
}

Integer& Integer::div(Integer& res, const Integer& n1, const long n2)
{
	if (isZero(n1)) return res = Integer::zero;
	//  if (isZero(n2)) {
	//    GivMathDivZero("[Integer::/]: division by zero");
	//  }
	int sgn = GMP__SGN(n2);
	mpz_tdiv_q_ui( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, GMP__ABS(n2));
	if (sgn <0) return res = -res;

	// if (n2>0)
	// mpz_fdiv_q_ui( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&n1.gmp_rep, n2);
	// else
	// mpz_cdiv_q_ui( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&n1.gmp_rep, n2);


	return res;
}


Integer& Integer::div(Integer& res, const Integer& n1, const int n2) {
	return div(res,n1,long(n2));
}

Integer& Integer::div(Integer& res, const Integer& n1, const unsigned long n2)
{
	if (isZero(n1)) return res = Integer::zero;
	//  if (isZero(n2)) {
	//    GivMathDivZero("[Integer::/]: division by zero");
	//  }
	mpz_tdiv_q_ui( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, n2);
	return res;
}

Integer& Integer::divexact  (Integer& q, const Integer& n1, const Integer& n2)
{
	if (isZero(n1)) return q = Integer::zero;
	//  if (isZero(n2)) {
	//    GivMathDivZero("[Integer::/]: division by zero");
	//  }
	mpz_divexact( (mpz_ptr)&(q.gmp_rep),
		      (mpz_ptr)&(n1.gmp_rep), (mpz_ptr)&(n2.gmp_rep)) ;
	return q;
}

Integer  Integer::divexact  (const Integer& n1, const Integer& n2)
{
	if (isZero(n1)) return Integer::zero;
	//  if (isZero(n2)) {
	//    GivMathDivZero("[Integer::/]: division by zero");
	//  }
	Integer q;
	mpz_divexact( (mpz_ptr)&(q.gmp_rep),
		      (mpz_ptr)&(n1.gmp_rep), (mpz_ptr)&(n2.gmp_rep)) ;
	return q;
}


Integer& Integer::operator /= (const Integer& n)
{
	//  if (isZero(n)) {
	//    GivMathDivZero("[Integer::/]: division by zero");
	//  }
	if (isZero(*this)) return *this;
	mpz_tdiv_q( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, (mpz_ptr)&n.gmp_rep) ;
	return *this;
}

Integer& Integer::operator /= (const unsigned long l)
{
	//  if (l ==0) {
	//    GivMathDivZero("[Integer::/]: division by zero");
	//  }
	if (isZero(*this)) return *this;
	mpz_tdiv_q_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, l);
	return *this;
}

Integer& Integer::operator /= (const long l)
{
	//  if (l ==0) {
	//    GivMathDivZero("[Integer::/]: division by zero");
	//  }
	if (isZero(*this)) return *this;
	int sgn = GMP__SGN(l);
	mpz_tdiv_q_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, GMP__ABS(l));
	if (sgn <0) mpz_neg( (mpz_ptr)&gmp_rep, (mpz_ptr)&(gmp_rep));
	return *this;
}


Integer Integer::operator / (const Integer& n) const
{
	//  if (isZero(n)) {
	//    GivMathDivZero("[Integer::/]: division by zero");
	//  }
	if (isZero(*this)) return Integer::zero;
	Integer res;
	mpz_tdiv_q( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, (mpz_ptr)&n.gmp_rep) ;

	// if (n>0)
		// mpz_fdiv_q( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, (mpz_ptr)&n.gmp_rep) ;
	// else
		// mpz_cdiv_q( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, (mpz_ptr)&n.gmp_rep) ;

	return res;
}

Integer Integer::operator / (const unsigned long l) const
{
	//  if (l ==0) {
	//    GivMathDivZero("[Integer::/]: division by zero");
	//  }
	if (isZero(*this)) return Integer::zero;
	Integer res;
	mpz_tdiv_q_ui( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, l);

	return res;
}

Integer Integer::operator / (const long l) const
{
	//  if (l ==0) {
	//    GivMathDivZero("[Integer::/]: division by zero");
	//  }
	if (isZero(*this)) return Integer::zero;
	Integer res;
	int sgn = GMP__SGN(l);
	mpz_tdiv_q_ui( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, GMP__ABS(l));
	if (sgn <0) return negin(res);
	// if (l>0)
		// mpz_fdiv_q_ui( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, l) ;
	// else
		// mpz_fdiv_q_ui( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, -l ) ;


	return res;
}

// -- Euclidian division
Integer& Integer::divmod(Integer& q, Integer& r, const Integer &a, const Integer &b)
{
	//  if (isZero(b)) {
	//    GivMathDivZero("[Integer::divide]: division by zero");
	//  }

	mpz_tdiv_qr( (mpz_ptr)&(q.gmp_rep), (mpz_ptr)&(r.gmp_rep),
		     (mpz_ptr)&(a.gmp_rep), (mpz_ptr)&(b.gmp_rep));
	// if (a>0)
		// mpz_fdiv_qr_ui( (mpz_ptr)&(q.gmp_rep),r
				// (mpz_ptr)&(a.gmp_rep), b);
	// else
		// mpz_cdiv_qr_ui( (mpz_ptr)&(q.gmp_rep),r
				// (mpz_ptr)&(a.gmp_rep), b);

	if (a<0 && r) {
		subin(q,(long)1) ;
		r = b - r;
	}


	return q;
}

Integer& Integer::divmod(Integer& q, long& r, const Integer& a, const long b)
{
	//  if (isZero(b)) {
	//    GivMathDivZero("[Integer::divide]: division by zero");
	//  }
	// int sgn = GMP__SGN(b);
	r = mpz_tdiv_q_ui( (mpz_ptr)&(q.gmp_rep),
			   (mpz_ptr)&(a.gmp_rep), GMP__ABS(b));
	// if (sgn <0) return negin(q);
	// if (a>0)
		// mpz_fdiv_qr_ui( (mpz_ptr)&(q.gmp_rep),r
				// (mpz_ptr)&(a.gmp_rep), b);
	// else
		// mpz_cdiv_qr_ui( (mpz_ptr)&(q.gmp_rep),r
				// (mpz_ptr)&(a.gmp_rep), b);
	if (a<0 && r) {
		subin(q,(long)1) ;
		r = b - r ;
	}



	return q;
}

Integer& Integer::divmod(Integer& q, unsigned long& r, const Integer& a, const unsigned long b)
{
	//  if (isZero(b)) {
	//    GivMathDivZero("[Integer::divide]: division by zero");
	//  }
	r = mpz_tdiv_q_ui( (mpz_ptr)&(q.gmp_rep),
			   (mpz_ptr)&(a.gmp_rep), b);

	if (a<0 && r) {
		subin(q,(long)1) ;
		r = b - r;
	}

	return q;
}

}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen
