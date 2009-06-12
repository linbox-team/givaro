#include <iostream>
#include <stdlib.h>
#include <givaro/givgfq.h>
#include <givaro/givpoly1factor.h>
#include <givaro/givtimer.h>

using namespace std;


int main(int argc, char** argv)
{
  GFqDom<long>::Residu_t MOD;
  if (argc > 1) 	
	  MOD = atoi(argv[1]);
  else
	  std::cin >> MOD;
  unsigned long expo = 1;
  if (argc > 2) expo = atoi(argv[2]);

  GFqDom<long> F(MOD, expo);
  
  Poly1FactorDom<GFqDom<long>, Dense> FD(F,Indeter("X"));
  typedef Poly1FactorDom<GFqDom<long>, Dense>::Element Polys ;
  Polys P;
  FD.read( cin, P );
  std::vector<Polys> Lf;
  std::vector<unsigned long> Le;
  
        Timer tim; tim.clear(); tim.start();
        FD.CZfactor(Lf, Le, P);
        tim.stop();

  FD.write( cout, P ) << " is 1";
  std::vector<unsigned long>::const_iterator e = Le.begin();
  for(std::vector<Polys>::const_iterator i = Lf.begin(); i != Lf.end(); ++i, ++e) {
      FD.write(cout << " * (", *i) << ")";
      if (*e > 1) cout << "^" << *e;
  }
  std::cout << std::endl;
  std::cerr << tim << std::endl;

  return 0;
}












//         bool f;

//     Poly1FactorDom<GFqDom<long>, Dense>::Element W,D; 
//     FD.gcd(W,FD.diff(D,P),P);
//     Degree d, dP;
//     if (FD.degree(d,W) > 0) return 0;
//         // Distinct degree free ?
//     Poly1FactorDom<GFqDom<long>, Dense>::Element  Unit, G1; 
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



