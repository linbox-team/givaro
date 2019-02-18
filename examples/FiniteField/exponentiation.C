// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

/*! @file examples/FiniteField/exponentiation.C
 * @ingroup examples
 * @ingroup finitefields
 * @example examples/FiniteField/exponentiation.C
 * @brief NO DOC
 */

#include <iostream>
#include <givaro/givpower.h>
#include <givaro/modular.h>
#include <givaro/gfq.h>

using namespace Givaro;



int main(int argc, char ** argv) {

    {
        Modular<int32_t> Z13(13);   // modulo 13 over 32 bits
        Modular<int32_t>::Element a, c;
        Z13.init(a, 7);

        long l = 29;
        dom_power(c, a, l, Z13); // c = 7^29 modulo 13 by squaring

        std::cerr << "Within ";
        Z13.write( std::cerr );
        std::cerr << " : " << std::flush;

        // Separate output writing
        Z13.write( std::cout, a) << " ^ " << l << " = " << std::flush;
        Z13.write( std::cerr, c) << std::endl;
    }


    {
        typedef GFqDom<int>::Residu_t TT;
        int Mod = 13; int exponent = 1;
        GFqDom<int> GF13( (TT) Mod, (TT) exponent );  // finite field with 13 elements
        GFqDom<int>::Element a, c;

        GF13.init(a, 7);    // 7 modulo   13

        long l = 29;
        dom_power(c, a, l, GF13); // c = 7^29 modulo 13 by squaring

        // Writing all outputs in a single command line
        GF13.write( std::cerr << "Within " ) << " : " << std::flush;
        GF13.write( GF13.write(
                               std::cout, a) << " ^ " << l << " = ", c) << std::endl;

    }

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
