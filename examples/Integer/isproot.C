#include <iostream>
#include <givaro/givintnumtheo.h>
#include <givaro/givtimer.h>



int main(int argc, char** argv)
{
  IntNumTheoDom<> IP;
  IntNumTheoDom<>::element a,b;
  if (argc > 1) a = Integer(argv[1]); else std::cin >> a;
  if (argc > 2) b = Integer(argv[2]); else std::cin >> b;
  
        Timer tim; tim.clear(); tim.start();
        bool f = IP.is_prim_root(a,b);
        tim.stop();
	std::cout << f << std::endl;
	std::cerr << tim << std::endl;

  return 0;
}

