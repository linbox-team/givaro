#include <iostream>
#include "givaro/givintrsa.h"
#include "givaro/givrandom.h"
#include "givaro/givtimer.h"


int main(int argc, char** argv)
{
    Timer tim;
    tim.clear();
    
    // m will be p times q
    // Sizes are in number of int to compose the big integer
    long keysize;
    if (argc > 1)
        keysize = atoi( argv[1] );
    else 
        std::cin >> keysize;


    GivRandom generator;


    tim.start();
    IntRSADom<GivRandom> IR(keysize,true,generator);
    tim.stop();
    
    std::cout << IR.getn() << " " << IR.gete() << " " << IR.getd()  << std::endl;
    std::cerr << tim << std::endl;
   
    return 0;
}

