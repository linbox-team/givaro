// ==========================================================================
// Copyright(c)'1994-2010 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// file: gfqkronecker.h
// Time-stamp: <12 Apr 10 16:46:17 Jean-Guillaume.Dumas@imag.fr>
// ==========================================================================
//
/*! @file examples/FiniteField/kronecker.C
 * @ingroup examples
 * @ingroup finitefields
 * @example examples/FiniteField/kronecker.C
 * @brief NO DOC
 */
#include <iostream>
#include <givaro/givrandom.h>
#include <givaro/modular-integer.h>
#include <givaro/gfqkronecker.h>

using namespace Givaro;



int main(int argc, char ** argv) {

    unsigned long charac = (argc>1?atoi(argv[1]):3);
    unsigned long expo   = (argc>2?atoi(argv[2]):8);


    GFqKronecker<long,Integer>  GF(charac, expo);
    GFqKronecker<long,Integer>::Element a, b;
    GF.write(std::cerr << "After field init of ") << std::endl;

    GF.init(a, 7U);
    GF.write( std::cerr << "7 --> ", a ) << std::endl;

    GivRandom generator;

    GF.random(generator,b);

    GF.write( std::cerr << "b: ", b ) << std::endl;

    Modular<int32_t> Zp(charac);
    Poly1PadicDom<Modular<int32_t> > PAD(Zp ,"Y");
    Poly1PadicDom<Modular<int32_t> >::Element pol;
    PAD.radixdirect(pol, GF.convert(b), expo);
    PAD.write(std::cerr<< "b(Y): ", pol) << std::endl;


    Integer r;
    unsigned long shift = ceil(log( (double)charac)/log(2.));
    for(size_t i=0; i<5; ++i, ++shift) {
        GF.setShift( shift );

        GF.convert(r, b);
        std::cerr << "b kron(" << shift << "): " << r << std::endl;

        Modular<Integer> Zm( 1<<shift );
        Poly1PadicDom<Modular<Integer> > PmAD(Zm ,"Z");
        Poly1PadicDom<Modular<Integer> >::Element pol;
        PmAD.radixdirect(pol, r, expo);
        PmAD.write(std::cerr<< "b(" << (1<<shift) << "): ", pol) << std::endl;

    }

    --shift;
    /// Test arithmetic
    Modular<Integer> Zm( 1<<shift );
    Zm.write(std::cerr << "with shift: ") << std::endl;
    Poly1PadicDom<Modular<Integer> > PmAD(Zm ,"B");
    Poly1PadicDom<Modular<Integer> >::Element Ipol;


    GFqKronecker<long,Integer>::Element c,d,e,f;
    GF.random(generator,c);
    GF.random(generator,d);
    PAD.radixdirect(pol, GF.convert(c), expo);
    PmAD.radixdirect(Ipol, GF.convert(r,c), expo);
    PmAD.write(std::cerr<< "c:=", Ipol) << ';' << std::endl;
    PmAD.radixdirect(Ipol, GF.convert(r,d), expo);
    PmAD.write(std::cerr<< "d:=", Ipol) << ';' << std::endl;
    GF.random(generator,e);
    GF.random(generator,f);
    PmAD.radixdirect(Ipol, GF.convert(r,e), expo);
    PmAD.write(std::cerr<< "e:=", Ipol) << ';' << std::endl;
    PmAD.radixdirect(Ipol, GF.convert(r,f), expo);
    PmAD.write(std::cerr<< "f:=", Ipol) << ';' << std::endl;


    GFqKronecker<long,Integer>::Element dot;
    GF.mul(dot,c,d);
    GF.axpyin(dot,e,f);


    GF.write( std::cerr << "dot: ", dot ) << std::endl;
    PmAD.radixdirect(Ipol, GF.convert(r,dot), expo);
    PmAD.write(std::cerr<< "dot:=", Ipol) << ';' << std::endl;


    Integer Ic,Id,Ie,If;
    GF.convert(Ic, c);
    GF.convert(Id, d);
    GF.convert(Ie, e);
    GF.convert(If, f);


    Integer Idot = Ic*Id+Ie*If;

    PmAD.radixdirect(Ipol, Idot, 2*(expo));
    PmAD.write(std::cerr<< "Idot:=", Ipol) << ';' << std::endl;

    std::cerr << "Idot: " << Idot << std::endl;


    GFqKronecker<long,Integer>::Element res;
    GF.init(res, Idot);

    GF.write( std::cerr << "res: ", res ) << std::endl;



    if (! GF.areEqual(res, dot)) {
        std::cerr << "ERROR: incoherency" << std::endl;
    } else
        std::cerr << "passed." << std::endl;

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
