// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

/*! @file examples/Integer/primitiveelement.C
 * @ingroup examples
 * @ingroup integers
 * @example examples/Integer/primitiveelement.C
 * @brief NO DOC
 */
#include <iostream>
#include <givaro/givintnumtheo.h>
#include <givaro/givtimer.h>


using namespace Givaro;




int main(int argc, char** argv)
{
  IntNumTheoDom<> IP;
  IntNumTheoDom<>::Element a,pr;
  if (argc > 1) a = IntNumTheoDom<>::Element(argv[1]); else std::cin >> a;

        Timer tim; tim.clear(); tim.start();
	IP.prim_elem(pr, a);
        tim.stop();
	IntegerDom().write( std::cout, pr ) << std::endl;
	std::cerr << tim << std::endl;

  return 0;
}

