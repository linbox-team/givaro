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
    
    IntRSADom<GivRandom>::Rep n,e,d;
    if (argc > 3)
        n = Integer( argv[3] );
    else 
        std::cin >> e;
    if (argc > 4)
        e = Integer( argv[4] );
    else 
        std::cin >> e;
    if (argc > 5)
        d = Integer( argv[5] );
    else 
        std::cin >> d;
    
    IntRSADom<GivRandom> IR(n,e,d);

    std::ifstream TXT(argv[1]);
    if (!TXT) { std::cerr << "Error opening input file: " << argv[1] << std::endl; return -1; }
    std::ofstream OUT(argv[2]);
    if (!OUT) { std::cerr << "Error opening output file: " << argv[2] << std::endl; return -1; }

    
    tim.start();
    IR.decipher( OUT, TXT );
    OUT.close();
    TXT.close();
    tim.stop();

    std::cerr << tim << std::endl;
    
   
    return 0;
}

