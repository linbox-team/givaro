#include <iostream>
#include <stdlib.h>
#include <givaro/givinteger.h>
#include <givaro/givintprime.h>
#include <givaro/givtimer.h>


int main(int argc, char** argv)
{
    Integer m, p;
    if (argc > 1) m = Integer(argv[1]);
    IntPrimeDom IP;
    
    {
        Timer tim; tim.clear(); tim.start();
        int a = isperfectpower(m);
        tim.stop();
        std::cout << a << std::endl;
        std::cerr << tim << std::endl;
    }
    {
        Timer tim; tim.clear(); tim.start();
        int a = IP.isprimepower(p, m);
        tim.stop();
        if (a) std::cout << "is " << p << "^" << a << std::endl;
        else   std::cout << "not a prime power" << std::endl;

        std::cerr << tim << std::endl;
    }
    return 0;
}
