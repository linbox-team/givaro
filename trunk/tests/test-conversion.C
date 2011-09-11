// Copyright(c)'2010 by The Givaro group
// This file is part of Givaro.
// written by BB
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.


#include <givaro/givinteger.h>
using namespace Givaro;

int main()
{
	int qq = -1 ;
	Integer q(qq);
	int q1           = (int) q;
	long int q2      = (long int) q;
	long long int q3 = (long long int) q;
	if (q1 != qq)
		return 1 ;
	if (q2 != qq)
		return 2 ;
	if (q3 != qq)
		return 3 ;
	return 0 ;
}
