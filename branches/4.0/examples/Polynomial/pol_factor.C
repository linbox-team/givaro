// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
/*! @file examples/Polynomial/pol_factor.C
 * @ingroup examples
 * @ingroup polynomials
 * @example examples/Polynomial/pol_factor.C
 * @brief NO DOC
 */

#include <iostream>
#include <stdlib.h>
#include <givaro/gfq.h>
#include <givaro/givpoly1factor.h>
#include <givaro/givtimer.h>

using namespace std;

using namespace Givaro;



int main(int argc, char** argv)
{
  GFqDom<int64_t>::Residu_t MOD;
  if (argc > 1)
	  MOD = (GFqDom<int64_t>::Residu_t) atoi(argv[1]);
  else
	  std::cin >> MOD;
  uint64_t expo = 1;
  if (argc > 2) expo = (uint64_t)atoi(argv[2]);

  GFqDom<int64_t> F(MOD, expo);

  Poly1FactorDom<GFqDom<int64_t>, Dense> FD(F,Indeter("X"));
  typedef Poly1FactorDom<GFqDom<int64_t>, Dense>::Element Polys ;
  Polys P;
  FD.read( cin, P );
  std::vector<Polys> Lf;
  std::vector<uint64_t> Le;

        Timer tim; tim.clear(); tim.start();
        FD.CZfactor(Lf, Le, P);
        tim.stop();

  FD.write( cout, P ) << " is 1";
  std::vector<uint64_t>::const_iterator e = Le.begin();
  for(std::vector<Polys>::const_iterator i = Lf.begin(); i != Lf.end(); ++i, ++e) {
      FD.write(cout << " * (", *i) << ")";
      if (*e > 1) cout << "^" << *e;
  }
  std::cout << std::endl;
  std::cerr << tim << std::endl;

  return 0;
}




#if 0







		bool f;

	Poly1FactorDom<GFqDom<int64_t>, Dense>::Element W,D;
	FD.gcd(W,FD.diff(D,P),P);
	Degree d, dP;
	if (FD.degree(d,W) > 0) return 0;
		// Distinct degree free ?
	Poly1FactorDom<GFqDom<int64_t>, Dense>::Element  Unit, G1;
	FD.init(Unit, Degree(1), F.one);
	W.copy(Unit);
	FD.degree(dP,P); Degree dPo = (dP/2);

	f = 1;

	for(Degree dp = 1; dp <= dPo; ++dp) {
		FD.write(cout << "W: ", W) << endl ;
		FD.powmod(W, D.copy(W), MOD, P);
		FD.write(cout << "W^q: " , W) << endl ;
		FD.gcd (G1, FD.sub(D,W,Unit), P) ;
		FD.write(cout << "D: " , D) << endl ;
		FD.write(cout << "G1: " , G1) << endl ;
		if ( FD.degree(d,G1) > 0 ) { f = 0; break; }
	}


#endif

