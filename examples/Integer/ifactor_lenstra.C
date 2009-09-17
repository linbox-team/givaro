// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software. 
// see the COPYRIGHT file for more details.

#define GIVARO_LENSTRA
#include <iostream>
using namespace std;
#include <givaro/givinit.h>


#include <givaro/givintfactor.h>
#include <givaro/givtimer.h>

int main(int argc, char** argv)
{
    IntFactorDom<> IP;
#ifndef __GIVARO_GMP_NO_CXX
      IP.seeding();
      // std::cerr << "Seeding..." << std::endl;
#endif
    Integer m;
    if (argc > 1)
       m = Integer(argv[1]);
    else
        cin >> m;
    if (IP.islt(m,0) ) {
        cerr << "-";
        IP.negin(m);
   } 
    if (IP.islt(m,4))
        IP.write(cerr,m) << endl;
    else {
        Timer tim; tim.clear(); tim.start();
        IP.write(cerr,m) << endl;
        tim.stop();
        cerr << tim << endl;
    }
    return 0;
}

