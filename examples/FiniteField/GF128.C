// ========================================================== //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <12 Jun 15 18:38:55 Jean-Guillaume.Dumas@imag.fr>
// Thanks to Dieter Schuster
// ========================================================== //

/*! @file examples/FiniteField/GF128.C
 * @ingroup examples
 * @ingroup finitefields
 * @example examples/FiniteField/GF128.C
 * @brief NO DOC
 */
#include <givaro/gfq.h>
#include <givaro/givtimer.h>

using namespace Givaro;



int main(int argc, char** argv)
{
    GFqDom<int64_t> GF128(2, 7);
    GFqDom<int64_t>::Element b, c, gen;
    GF128.init(b, 5U);
    GF128.init(c, 3);
    GF128.write(std::cout, b) << std::endl;
    GF128.write(std::cout, c) << std::endl;

    GFqDom<int64_t>::Element f,g,h,j;
    GFqDom<int64_t> F2(2);
    Poly1Dom< GFqDom<int64_t>, Dense> Pol2(F2);
    Poly1Dom< GFqDom<int64_t>, Dense>::Element P, Q, R;
    Pol2.init(P,Degree(1));
    F2.init(P[0],1);
    F2.init(P[1],1);
    GF128.init(f, P);
    GF128.write(std::cout << "2-adic representation of 1+X is: ", f)
    << std::endl
    << " ... while its internal representation is: "
    << f << std::endl;

    GF128.write(std::cout << "Indeed, we are in ") <<std::endl;

    GF128.generator(gen);
    GF128.write(std::cout <<
                "In this field, the generator used is (in 2-adic): ", gen)
    << std::endl
    << " whose internal representation is "
    << gen << std::endl;


    Poly1PadicDom< GFqDom<int64_t>, Dense > Padic2(Pol2);
    //

    std::cout << "Irreducible (in 2-adic): "
    << GF128.irreducible() << std::endl;

    GF128.init(g, Padic2.radix( Q, Integer(5) ));
    GF128.write(std::cout
                << "2-adic representation of 1+X^2 is: ", g)
    << std::endl;

    GF128.init(h);
    GF128.add(h, g, f);
    GF128.write(std::cout
                << "2-adic representation of X+X^2 is: ", h)
    << std::endl;

    GF128.mul(h, g, f);
    GF128.write(std::cout
                << "2-adic representation of 1+X+X^2+X^3 is: ", h)
    << std::endl;

    GF128.div(h, g, f);
    GF128.write(std::cout
                << "2-adic representation of 1+X is: ", h)
    << std::endl;

    GF128.init(j, Padic2.radix( Q, Integer(213) ));
    GF128.write(std::cout
                << "2-adic representation of the moding out of X^7+X^6+X^4+X^2+1 by the irreducible is: ", j)
    << std::endl;

    return 0;

}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
