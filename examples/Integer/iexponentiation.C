// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

/*! @file examples/Integer/iexponentiation.C
 * @ingroup examples
 * @ingroup integers
 * @example examples/Integer/iexponentiation.C
 * @brief NO DOC
 */
#include <iostream>
#include <givaro/modular-integer.h>
#include <givaro/givpower.h>
#include <givaro/givtimer.h>

using namespace Givaro;


// Modular exponentiation
// argv[1] : a
// argv[2] : e
// argv[3] : p
// res ---> a^e % p
int main(int argc, char ** argv) {

    {
        Integer p(argv[3]);
        Modular<Integer> Zp( p );
        Modular<Integer>::Element a, b;
        unsigned long e = (unsigned long)atoi(argv[2]) ;
        Zp.init(a, Integer(argv[1]));
        Zp.init(b);

        Timer tim;tim.clear();tim.start();
        dom_power(b, a, (long)e, Zp);
        tim.stop();

        Zp.write( std::cout << '(', a) << '^' << e << ") % " << p << '=' << std::flush;
        Zp.write( std::cerr, b) << std::endl;

        std::cerr << tim << std::endl;

    }




    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
