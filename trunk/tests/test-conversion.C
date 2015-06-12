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

	int64_t q2 = (int64_t) q;
	if (q2 != qq)
		return err = 2;

	int64_t q3 = (int64_t) q;
	if (q3 != qq)
		return err = 3 ;

	/*  test unsigned versions */
	/*  cast towards unsigned consider only the absolute value */
	int64_t uqq = (int64_t) qq;
	int64_t q4  = (int64_t) (int) q;
// 	std::cout << q4 << std::endl;
// 	std::cout << uqq << std::endl;
	if (q4 != uqq)
		return err = 4;

	uint64_t lqq = (uint64_t) qq;
//         std::cerr << "q  : " << q << std::endl;
	uint64_t q5 = (uint64_t) (int64_t) q;
//         std::cerr << "lqq: " << lqq << std::endl;
//         std::cerr << "q5 : " << q5 << std::endl;
	if (q5 != lqq)
		return err = 5;

	uint64_t luqq = (uint64_t) qq;
	uint64_t q6 = (uint64_t) (int64_t) q;
	if (q6 != luqq)
		return err = 6 ;
#


	/*  test unsigned versions */
	/*  cast towards unsigned consider only the absolute value */
	int64_t vqq = (int64_t) qq;
	int64_t q7  = (int64_t) q;
// 	std::cout << q7 << std::endl;
// 	std::cout << uqq << std::endl;
	if (q7 != -vqq)
		return err = 7;

	uint64_t q8 = (uint64_t) q;
	if (q8 != (uint64_t)(-vqq))
		return err = 8;

	uint64_t q9 = (uint64_t) q;
	if (q9 != (uint64_t)(-vqq))
		return err = 9 ;

	return err ;

}

int main()
{
	int err = 0;
	err  = testBasicConversion();
	return err ;

}
