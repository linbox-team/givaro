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

#if 0
int test_int_cstor()
{

	int a      = 7;
	Rationel A(a);
	if (A != a)
		return false ;
}
#endif

int main()
{
	int a      = 7;
	long int b = -6 ;
	unsigned long int c = 100 ;
	Rationel A(a);
	Rationel B(3);
	Rationel C(A);
	Rationel D(b);
	Rationel E(-5L);
	Rationel F(c);
	Rationel G(10UL);
	Rationel H(0.3);
	Rationel I(1,2);
	Rationel J(7UL,-2L);
	std::cout << "xxx" << std::endl;
	Rationel K(a,b);
	std::cout << "xxx" << std::endl;
	Rationel L(a,1);
	std::istringstream s ;
	s.str("4/5");
	Rationel M ;
	s >> M ;
	Rationel N(b,-15);
	std::cout << "---------------------------" << std::endl;
	std::cout << a       << '=' << A          << std::endl;
	std::cout << "3 = "  << B   << std::endl;
	std::cout << a       << "=" << C          << std::endl;
	std::cout << b       << "=" << D          << std::endl;
	std::cout << "-5 = " << E   << std::endl;
	std::cout << "---------------------------" << std::endl;
	std::cout << c       << "=" << F          << std::endl;
	std::cout << "10="   << G   << std::endl;
	std::cout << "0.3="  << H   << "="        << H.reduce() << std::endl;
	std::cout << "1/2="  << I   << std::endl;
	std::cout << "-7/2=" << J   << std::endl;
	std::cout << "---------------------------" << std::endl;
	std::cout << a       << "/" << b          << "="        << K          << std::endl;
	K.reduce();
	std::cout << a       << "/" << b          << "="        << K          << std::endl;
	std::cout << a       << "=" << L          << std::endl;
	std::cout << "4/5="  << M   << std::endl;
	std::cout << J       << "+" << A          << "=" ;
	Rationel::addin(J,A);
	std::cout << J << std::endl;
	std::cout << b       << "/-15="        << N          << std::endl;
	N.reduce();
	std::cout << b       << "/-15="        << N          << std::endl;
	s.str("4/-5");
	s >> M ;
	std::cout << "4/-5="  << M   << std::endl;


	return 0 ;
}
