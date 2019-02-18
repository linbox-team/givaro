// =================================================================== //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <21 Feb 12 17:02:32 Jean-Guillaume.Dumas@imag.fr>
// Givaro : Modular square roots
// =================================================================== //

/*! @file examples/Integer/ModularSquareRoot.C
 * @ingroup examples
 * @ingroup integers
 * @example examples/Integer/ModularSquareRoot.C
 * @brief NO DOC
 */
#include <iostream>
#include <stdlib.h>
#include <givaro/givintsqrootmod.h>
#include <givaro/givtimer.h>


#include <givaro/givpoly1factor.h>
#include <givaro/modular-integer.h>

using namespace Givaro;

// Algorithm 3.34 (Square Root Mod p) of
// Handbook of Applied Cryptography
// by Menezes, van Oorschot, Vanstone

int main(int argc, char** argv)
{
    Integer a(argv[1]), n(argv[2]);
    IntSqrtModDom<> ISM;

    {

        std::cerr << "a: " << a << std::endl;
        std::cerr << "n: " << n << std::endl;

        Integer::seeding ( (unsigned long)BaseTimer::seed ());

        Integer r;
        Timer chrono; chrono.start();
        ISM.sqrootmod(r,a,n);
        chrono.stop();
        std::cout << r << std::endl;
        std::cerr << chrono << std::endl;

        std::cerr << "Check, (" << r << ")^2 mod " << n << " = " << ( (r*r)%n) << std::endl;
    }

    if (ISM.isprime(n)) {
        std::cout << "Using polynomial factorization : " << std::endl;
        typedef Modular<Integer> Field;
        typedef Poly1FactorDom<Field,Dense> Polys;
        typedef Polys::Element Polynomial;
        Field F(n); Polys Pol(F, "X");

        Polynomial quad, root;
        Pol.init(quad, Degree(2));
        F.init(quad[0],a); F.negin(quad[0]);

        Timer chrono; chrono.start();
        if (Pol.is_irreducible(quad)) {
            std::cerr << a << " is not a quadratic residue mod " << n << std::endl;
            chrono.stop();
            std::cerr << chrono << std::endl;
            return 0;
        }

        Pol.SplitFactor(root, quad, Degree(1));
        chrono.stop();

        Pol.divin(root, root[1]);


        Pol.write(std::cout, root) << std::endl;
        std::cerr << chrono << std::endl;

        std::cerr << "Check, (" << root[0] << ")^2 mod " << n << " = " << ( (root[0]*root[0])%n) << std::endl;
    }


    return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
