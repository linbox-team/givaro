// Copyright(c)'2010 by The Givaro group
// This file is part of Givaro.
// written by BB
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.


#include <givaro/givinteger.h>
#include <math.h>

using namespace Givaro;

int testBasicConversion()
{
	int err = 0;
	/* test small neg int */
	int qq = -1 ;
	Integer q(qq);
	int q1  = (int) q;
	if (q1 != qq)
		err  = 1;

	long int q2 = (long int) q;
	if (q2 != qq)
		err  = 2;

#ifndef __GIVARO__DONOTUSE_longlong__
	long long int q3 = (long long int) q;
	if (q3 != qq)
		err = 3 ;
#endif

	/*  test unsigned versions */
	unsigned int uqq = (unsigned int) qq;
	unsigned int q4  = (unsigned int) q;
	std::cout << q4 << std::endl;
	std::cout << uqq << std::endl;
	if (q4 != uqq)
		err = 4;

	unsigned long int q5 = (unsigned long int) q;
	if (q5 != uqq)
		err = 5;

#ifndef __GIVARO__DONOTUSE_longlong__
	unsigned long long int q6 = (unsigned long long int) q;
	if (q6 != uqq)
		err = 6 ;
#endif


	return (err != 0) ;

}

int main()
{
	int err = 0;
	err  = testBasicConversion();
	return err ;

}
