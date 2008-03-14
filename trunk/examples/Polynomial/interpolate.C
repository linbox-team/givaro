#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <givaro/givzpz.h>
#include <givaro/givpoly1factor.h>
#include <givaro/givtimer.h>
#include <givaro/givinterp.h>

using namespace std;


int main(int argc, char** argv)
{
  ZpzDom<Std32>::Residu_t MOD;
  if (argc > 2) 	
	  MOD = atoi(argv[2]);
  else
	  std::cin >> MOD;

  ZpzDom<Std32> F(MOD);

  
  Interpolation< ZpzDom<Std32> > FD(F,Indeter("X"));
  Interpolation< ZpzDom<Std32> >::Element nouv, prec;
  int EarlyTerm = 0, Bound = 5;

  ZpzDom<Std32>::Element x, f;

  std::ifstream input (argv[1]);
  F.read(input, x);
  F.read(input, f);
  FD(x,f);
  prec = FD.interpolator();

  Timer tim; tim.clear(); tim.start();
  while(! F.read(input, x).eof() ) {
      F.read(input, f);
      FD(x,f);
      nouv = FD.interpolator();
      if (FD.areEqual(nouv, prec)) {
          if (++EarlyTerm > Bound) { std::cerr << "EarlyTerminated" << std::endl ; break; }
      } else
          EarlyTerm = 0;
      prec = nouv;
  }
  tim.stop();

  FD.write( std::cout, FD.interpolator()  ) << std::endl;
  std::cerr << tim << std::endl;

  return 0;
}




