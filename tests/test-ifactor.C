// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software. 
// written by BB
// see the COPYRIGHT file for more details.

#include <iostream>
using namespace std;
#include <givaro/givinit.h>
#include <givaro/givintfactor.h>
#include <givaro/givintprime.h>

#define _ITERS 100
//#define LOOPS 0

int test(const IntFactorDom<> & IP, const Integer & m)
{/*{{{*/

	Integer f  = 1 ;
	Integer r = 1 ;

	//IP.factor(f,m,LOOPS);
	IP.factor(f,m) ; // ne teste que Lenstra ou Pollard selon que que GIVARO_LENSTRA est définie ou non

	Integer::mod(r,m,f);
	if (r || !probab_prime(f) )  return -1 ;


	return 0 ;
}/*}}}*/

int main(int argc, char** argv)
{/*{{{*/
	IntFactorDom<> IP;
	Integer m;
	int err = 0 ;
	long int a = 3;
	long int b = 50 ;

	for (size_t i = 0 ; i < _ITERS ; ++i)
	{/*{{{*/
		if (!(i%25)) std::cout << i << "..." ;
		m = Integer::random_between(a,b);
		err = test(IP,m);
		if (err) break ;

	}/*}}}*/
	std::cout << _ITERS << std::endl;

	if (err) return err ;

	Integer p,q ;
	a = 29 ;
	b = 30 ;
	for (size_t i = 0 ; i < _ITERS ; ++i)
	{/*{{{*/
		if (!(i%25)) std::cout << i << "..." ;

		p = Integer::random_between(a,b);
		IP.nextprimein(p);
		q = Integer::random_between(a,b);
		IP.nextprimein(q);

		err = test(IP,m);
		if (err) break ;

		m = p*q ;

	}/*}}}*/
	std::cout << _ITERS << std::endl;

	if (err) return err ;

	return 0;
}/*}}}*/

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen

