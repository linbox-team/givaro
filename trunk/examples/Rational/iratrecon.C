#include <iostream>
#include <stdlib.h>
#include <givaro/givrational.h>
#include <givaro/givtimer.h>
#include <givaro/givinit.h>         // Givaro initialization



int main(int argc, char** argv)
{
    Integer f,m,k;
    if (argc > 1) f = Integer(argv[1]); else std::cin >> f;
    if (argc > 2) m = Integer(argv[2]); else std::cin >> m;
  
    RationalDom RD;
    Rational rec; 
  
    Timer tim; tim.clear(); tim.start();
    if (argc > 3) 
        RD.ratrecon(rec,f,m, Integer(argv[3]) );
    else 
        RD.ratrecon(rec,f,m);
    tim.stop();
    std::cout << rec.nume() << "/" << rec.deno() << " = " << f << " mod " << m << std::endl;
    std::cerr << tim << std::endl;
  
    return 0;
}

