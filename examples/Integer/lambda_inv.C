#include <iostream>
#include <givaro/givintnumtheo.h>
#include <givaro/givtimer.h>

// Lambda function : order of a primitive invertible
//                (invertible element of maximal orbit size)
//


int main(int argc, char** argv)
{
  IntNumTheoDom<> IP;
  IntNumTheoDom<>::element a,pr;
  if (argc > 1) a = IntNumTheoDom<>::element(argv[1]); else std::cin >> a;
  
        Timer tim; tim.clear(); tim.start();
	IP.lambda_inv(pr, a);
        tim.stop();
	IntegerDom().write( std::cout, pr ) << std::endl;
	std::cerr << tim << std::endl;

  return 0;
}

