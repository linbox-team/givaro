#include <iostream>
#include <fstream>
#include "givaro/givintrsa.h"
#include "givaro/givrandom.h"
#include "givaro/givtimer.h"


int main(int argc, char** argv)
{
    Timer tim;
    tim.clear();
    
    IntRSADom<GivRandom>::Rep m,k;
    if (argc > 3)
        m = Integer( argv[3] );
    else 
        std::cin >> m;
    if (argc > 4)
        k = Integer( argv[4] );
    else 
        std::cin >> k;
    
    IntRSADom<GivRandom> IR(m,k);

    tim.start();
    std::ifstream TXT(argv[1]);
    std::ofstream OUT(argv[2]);
    IR.encipher( OUT, TXT );
    OUT.close();
    TXT.close();
    tim.stop();

    std::cerr << tim << std::endl;
    
   
    return 0;
}

