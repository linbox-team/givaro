// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software. 
// see the COPYRIGHT file for more details.

#include <iostream>
#include <stdlib.h>
#include <givaro/givgfq.h>
#include <givaro/givpoly1factor.h>
#include <givaro/givtimer.h>

int main(int argc, char** argv) {
  GFqDom<long>::Residu_t MOD;
  if (argc > 1) 	
	  MOD = atoi(argv[1]);
  else
	  std::cin >> MOD;
  unsigned long expo = 1;
  if (argc > 2) expo = atoi(argv[2]);

  GFqDom<long> F(MOD, expo);
  
  Poly1FactorDom<GFqDom<long>, Dense> FD(F,Indeter("X"));
  Poly1FactorDom<GFqDom<long>, Dense>::Element P, IXE;
  FD.init(IXE,Degree(1),F.one);
  FD.read( std::cin, P );

        Timer tim; tim.clear(); tim.start();
        bool f = FD.is_prim_root(IXE, P );
        tim.stop();

  F.write( FD.write( std::cout, P ) << " is " << (f?"":"not ") << "primitive in " ) << endl;

      // std::cout << f << std::endl;
  std::cerr << tim << std::endl;

  return 0;
}












//         bool f;

//     Poly1FactorDom<GFqDom<long>, Dense>::element W,D; 
//     FD.gcd(W,FD.diff(D,P),P);
//     Degree d, dP;
//     if (FD.degree(d,W) > 0) return 0;
//         // Distinct degree free ?
//     Poly1FactorDom<GFqDom<long>, Dense>::element  Unit, G1; 
//     FD.init(Unit, Degree(1), F.one);
//     W.copy(Unit);
//     FD.degree(dP,P); Degree dPo = (dP/2);

//     f = 1;

//     for(Degree dp = 1; dp <= dPo; ++dp) {
//         FD.write(cout << "W: ", W) << endl ;
//         FD.powmod(W, D.copy(W), MOD, P);
//         FD.write(cout << "W^q: " , W) << endl ;
//         FD.gcd (G1, FD.sub(D,W,Unit), P) ;
//         FD.write(cout << "D: " , D) << endl ;
//         FD.write(cout << "G1: " , G1) << endl ;
//         if ( FD.degree(d,G1) > 0 ) { f = 0; break; }
//     }




