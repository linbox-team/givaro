// Copyright(c)'2011 by The Givaro group
// This file is part of Givaro.
// written by BB
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.


#include <gmp++/gmp++_rat.h>
#include <iostream>
#include <sstream>

using namespace Givaro ;

int main()
{
	int a      = 1;
	long int b = 1 ;
	unsigned long int c = 1 ;
	Rationel A(a);
	Rationel B(1);
	Rationel C(A);
	Rationel D(b);
	Rationel E(1L);
	Rationel F(c);
	Rationel G(1UL);
	Rationel H(0.3);
	Rationel I(1,2);
	Rationel J(1UL,2L);
	Rationel K(a,b);
	Rationel L(a,1);
	std::istringstream s ;
	s.str("4/5");
	Rationel M ;
	s >> M ;
	std::cout << A << std::endl;
	std::cout << B << std::endl;
	std::cout << C << std::endl;
	std::cout << D << std::endl;
	std::cout << E << std::endl;
	std::cout << F << std::endl;
	std::cout << G << std::endl;
	std::cout << H << std::endl;
	std::cout << H.reduce() << std::endl;
	std::cout << I << std::endl;
	std::cout << J << std::endl;
	std::cout << K << std::endl;
	std::cout << L << std::endl;
	std::cout << M << std::endl;


	return 0 ;
}
