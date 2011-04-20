// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

/*! @file examples/Integer/ispower.C
 * @ingroup examples
 * @ingroup integers
 * @example examples/Integer/ispower.C
 * @brief NO DOC
 */
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
