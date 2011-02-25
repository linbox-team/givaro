// Copyright(c)'1994-2011 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// written by BB
// see the COPYRIGHT file for more details.

/*! @file tests/test-integer.h
 * @ingroup tests
 * tests integer.h fucntions not tested elsewhere.
 * @test tests integer.h fucntions not tested elsewhere.
 */

#include "gmp++/gmp++.h"

#define _GIV_REPET 10

int test_axpy_unit(const Integer & a, const Integer & b, const Integer & x)
{
	Integer y,z ;
	/*  axpy */
	z = a*x + b ;
	Integer::axpy(y,a,x,b);
	if (y != z)
		return -1 ;

	y = b ;
	Integer::axpyin(y,a,x);
	if (y != z)
		return -2 ;

	y = b ;
	Integer::axpy(y,a,x,y);
	if (y != z)
		return -3 ;

	/* axmy */
	z = a*x - b ;
	Integer::axmy(y,a,x,b);
	if (y != z)
		return -4 ;

	y = b ;
	Integer::axmyin(y,a,x);
	if (y != z)
		return -5 ;

	y = b ;
	Integer::axmy(y,a,x,y);
	if (y != z)
		return -6 ;

	/* maxpy */
	z =  b - a * x ;
	Integer::maxpy(y,a,x,b);
	if (y != z)
		return -7 ;

	y = b ;
	Integer::maxpyin(y,a,x);
	if (y != z)
		return -8 ;

	y = b ;
	Integer::maxpy(y,a,x,y);
	if (y != z)
		return -9 ;

	return 0 ;

}

int test_axpy()
{
	Integer a,x,b ;
	int res = 0 ;
	int repet = _GIV_REPET ;
	while (--repet) {
		/* axpy */
		a = Integer::random<false>();
		b = Integer::random<false>();
		x = Integer::random<false>();
		res = test_axpy_unit(a,b,x);
		if (res)
			return res ;

		a = 0 ;
		res = test_axpy_unit(a,b,x);
		if (res)
			return res ;

		a = Integer::random<false>();
		b = 0 ;
		res = test_axpy_unit(a,b,x);
		if (res)
			return res ;

		b = Integer::random<false>();
		x = 0 ;
		res = test_axpy_unit(a,b,x);
		if (res)
			return res ;
	}

	return res ;
}

int test_mul_unit(const Integer & a, const Integer & b)
{
	Integer y,z ;
	/*  axpy */
	z = a*b ;
	Integer::mul(y,a,b);
	if (y != z)
		return -1 ;

	y = b ;
	Integer::mulin(y,a);
	if (y != z)
		return -2 ;

	y = a ;
	Integer::mul(y,y,b);
	if (y != z)
		return -3 ;

	y = b ;
	Integer::mul(y,a,y);
	if (y != z)
		return -3 ;

	z = b * b ;
	y = b ;
	Integer::mulin(y,y);
	if (y != z)
		return -5 ;

	return 0 ;

}

int test_mul()
{
	Integer a,b ;
	int res = 0 ;
	int repet = _GIV_REPET ;
	while (--repet) {
		/* axpy */
		a = Integer::random<false>();
		b = Integer::random<false>();
		res = test_mul_unit(a,b);
		if (res)
			return res ;

		a = 0 ;
		res = test_mul_unit(a,b);
		if (res)
			return res ;

		a = Integer::random<false>();
		b = 0 ;
		res = test_mul_unit(a,b);
		if (res)
			return res ;

	}

	return res ;
}

//! @todo test gcd...

int main ()
{
	int res = 0 ;
	res = test_axpy();
	if (res)
		return res ;

	res = test_mul();
	if (res)
		return res ;


	return res ;
}
