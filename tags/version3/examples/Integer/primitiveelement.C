#include <iostream>
#include <givaro/givintnumtheo.h>
#include <givaro/givtimer.h>



int main(int argc, char** argv)
{
  IntNumTheoDom<> IP;
  IntNumTheoDom<>::element a,pr;
  if (argc > 1) a = IntNumTheoDom<>::element(argv[1]); else std::cin >> a;
  
        Timer tim; tim.clear(); tim.start();
	IP.prim_elem(pr, a);
        tim.stop();
	IntegerDom().write( std::cout, pr ) << std::endl;
	std::cerr << tim << std::endl;

  return 0;
}

