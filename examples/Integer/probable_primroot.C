#include <iostream>
#include <givaro/givintnumtheo.h>
#include <givaro/givtimer.h>

//  Polynomial-time generation of primitive roots               
//  L is number of loops of Pollard partial factorization of n-1
//  10,000,000 gives at least 1-2^{-40} probability of success
//  [Dubrois & Dumas, Industrial-strength primitive roots]
//  Returns the probable primitive root and the probability of error.

int main(int argc, char** argv)
{
  IntNumTheoDom<> IP;
  IP.seeding( BaseTimer::seed() );

  double error;
  IntNumTheoDom<>::element a,pr;
  if (argc > 1) a = IntNumTheoDom<>::element(argv[1]); else std::cin >> a;

      // Default is partial factorization up to factors of at least 12 digits.
  Timer tim; tim.clear(); 
  tim.start();
  IP.probable_prim_root(pr, error, a );
  tim.stop();
  
  IntegerDom().write( std::cout << "Prim root   : ", pr ) << ", correct with probability at least : 1-" << error << std::endl;
  std::cerr << tim << std::endl;


  std::cerr << "Now checking primitivity, this may take some time (complete factorization of n-1) ...";

#define GIVARO_LENSTRA

  tim.start();
  if ( IP.isorder(a-1, pr, a) ) {
      tim.stop();
      std::cerr << "... Pimitivity checked" << std::endl;
      std::cerr << tim << std::endl;
  }
  else {
      tim.stop();
      std::cerr << "... WARNING : FAILURE" << std::endl;
      std::cerr << tim << std::endl;
  }
     
      

  return 0;
}

