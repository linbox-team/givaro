// Copyright(c)'2011 by The Givaro group
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
	const int qq = -1 ;
	const Integer q(qq);
	int q1  = (int) q;
	if (q1 != qq)
		return err = 1;

	long int q2 = (long int) q;
	if (q2 != qq)
		return err = 2;

#ifndef __GIVARO__DONOTUSE_longlong__
	long long int q3 = (long long int) q;
	if (q3 != qq)
		return err = 3 ;
#endif

	/*  test unsigned versions */
	/*  cast towards unsigned consider only the absolute value */
	unsigned int uqq = (unsigned int) qq;
	unsigned int q4  = (int) q;
// 	std::cout << q4 << std::endl;
// 	std::cout << uqq << std::endl;
	if (q4 != uqq)
		return err = 4;

	unsigned long lqq = (unsigned long) qq;
//         std::cerr << "q  : " << q << std::endl;
	unsigned long int q5 = (long int) q;
//         std::cerr << "lqq: " << lqq << std::endl;
//         std::cerr << "q5 : " << q5 << std::endl;
	if (q5 != lqq)
		return err = 5;

#ifndef __GIVARO__DONOTUSE_longlong__
	unsigned long long luqq = (unsigned long long) qq;
	unsigned long long int q6 = (long long int) q;
	if (q6 != luqq)
		return err = 6 ;
#endif


	/*  test unsigned versions */
	/*  cast towards unsigned consider only the absolute value */
	unsigned int vqq = (unsigned int) qq;
	unsigned int q7  = (unsigned int) q;
// 	std::cout << q7 << std::endl;
// 	std::cout << uqq << std::endl;
	if (q7 != -vqq)
		return err = 7;

	unsigned long int q8 = (unsigned long int) q;
	if (q8 != -vqq)
		return err = 8;

#ifndef __GIVARO__DONOTUSE_longlong__
	unsigned long long int q9 = (unsigned long long int) q;
	if (q9 != -vqq)
		return err = 9 ;
#endif

	return err ;

}

int main()
{
	int err = 0;
	err  = testBasicConversion();
	return err ;

}
