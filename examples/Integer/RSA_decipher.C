// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software. 
// see the COPYRIGHT file for more details.

#include <iostream>
#include <fstream>
#include "givaro/givintrsa.h"
#include "givaro/givrandom.h"
#include "givaro/givtimer.h"

// RSA, in CBC mode, deciphering of files



int main(int argc, char** argv)
{
    Timer tim;
    tim.clear();
    
    IntRSADom<GivRandom>::Rep m,k,u;
    if (argc > 3)
        m = Integer( argv[3] );
    else 
        std::cin >> m;
    if (argc > 4)
        k = Integer( argv[4] );
    else 
        std::cin >> k;
    if (argc > 5)
        u = Integer( argv[5] );
    else 
        std::cin >> u;
    
    IntRSADom<GivRandom> IR(m,k,u);

    tim.start();
    std::ifstream TXT(argv[1]);
    std::ofstream OUT(argv[2]);
    IR.decipher( OUT, TXT );
    OUT.close();
    TXT.close();
    tim.stop();

    std::cerr << tim << std::endl;
    
   
    return 0;
}

