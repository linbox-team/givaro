// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

/*! @file examples/FiniteField/zpz_atomic.C
 * @ingroup examples
 * @ingroup finitefields
 * @example examples/FiniteField/zpz_atomic.C
 * @brief NO DOC
 */
#include <stdio.h>
#include <stdlib.h>
#include <givaro/givtimer.h>
#include <givaro/givrandom.h>
#include <givaro/givzpz.h>

using namespace Givaro;



#ifndef GIVMIN
#define GIVMIN(a,b) (((a)<(b))?(a):(b))
#endif

typedef ZpzDom<Std32> Domain;
typedef Domain::Element Modulo;

#ifndef NB
#define NB 1000000
#endif
#ifndef TAILLE
#define TAILLE 256
#endif

int main(int argc, char ** argv)
{

    GivRandom generator;
    long P (32749);     // argv[1] : Modulus
    // int TAILLE (256);   // argv[2] : Vector size
    // long NB (500000);    // argv[3] : Number of iterations
    int offset = 0;
    if (argc > ++offset) P = atoi( argv[ offset ] );
    // if (argc > ++offset) TAILLE = atoi( argv[ offset ] );
    // if (argc > ++offset) NB = atoi( argv[ offset ] );

    Timer inver;
    inver.clear();
    inver.start();
    Domain GFq((Domain::Element)P);

    int i;
    std::cout << "."<< std::flush;

    Modulo z1, z2, z23, z3;

//    long seuil = GIVMIN(P*2,TAILLE);

    std::cout << "."<< std::flush;

        GFq.random(generator,z1) ;
        GFq.nonzerorandom(generator,z2) ;
        GFq.random(generator,z23) ;

    std::cout << "."<< std::flush;
    inver.stop();

    std::cout << "Init et Tableaux des  modulo " << P << " :\n" << inver << std::endl << std::flush;

    double coef = (double)NB*TAILLE / 1e6;

    int k;
    Timer tim;tim.clear();tim.start();
    for (k=0; k<NB; ++k)  for(i=0; i<TAILLE; ++i)
        GFq.mul(z3, z1, z2);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " Mul: " << coef / tim.usertime()
         << "Mop/s\n" << tim << ", ex: " << z3 << std::endl << std::flush;

    tim.clear();tim.start();
    for (k=0; k<NB; ++k)  for(i=0; i<TAILLE; ++i)
        GFq.add(z3, z1, z2);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " Add: " << coef / tim.usertime()
         << "Mop/s\n" << tim << ", ex: " << z3 << std::endl << std::flush;

    tim.clear();tim.start();
    for (k=0; k<NB; ++k)  for(i=0; i<TAILLE; ++i)
        GFq.axpyin(z3, z1, z2);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " MulAdd IN place: " << 2*coef / tim.usertime()
         << "Mop/s\n" << tim << ", ex: " << z3 << std::endl << std::flush;

    tim.clear();tim.start();
    for (k=0; k<NB; ++k)  for(i=0; i<TAILLE; ++i)
        GFq.axpy(z3, z1, z2, z23);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " MulAdd: " << 2*coef / tim.usertime()
         << "Mop/s\n" << tim << ", ex: " << z3 << std::endl << std::flush;

    tim.clear();tim.start();
    for (k=0; k<NB; ++k)  for(i=0; i<TAILLE; ++i)
        GFq.sub(z3, z1, z2);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " Sub: " << coef / tim.usertime()
         << "Mop/s\n" << tim << ", ex: " << z3 << std::endl << std::flush;
#if 0
    tim.clear();tim.start();
    for (k=0; k<NB; ++k)  for(i=0; i<TAILLE; ++i)
        GFq.div(z3, z1, z2);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " Div: " << coef / tim.usertime()
         << "Mop/s\n" << tim << ", ex: " << z3 << std::endl << std::flush;
#endif

    tim.clear();tim.start();
    for (k=0; k<NB; ++k) for(i=0; i<TAILLE; ++i)
        GFq.neg(z3, z1);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " Neg: " << (coef / tim.usertime())
         << "Mop/s\n" << tim << ", ex: " << z3 << std::endl << std::flush;

            GFq.random(generator,z1) ;
            GFq.nonzerorandom(generator,z2) ;
            GFq.random(generator,z23) ;


    tim.clear();tim.start();
    for (k=0; k<NB; ++k) for(i=0; i<TAILLE; ++i)
        GFq.add(z3, z1, z2);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " Add: " << coef / tim.usertime()
         << "Mop/s\n" << tim << ", ex: " << z3 << std::endl << std::flush;

    tim.clear();tim.start();
    for (k=0; k<NB; ++k) for(i=0; i<TAILLE; ++i)
        GFq.sub(z3, z1, z2);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " Sub: " << coef / tim.usertime()
         << "Mop/s\n" << tim << ", ex: " << z3 << std::endl << std::flush;

    tim.clear();tim.start();
    for (k=0; k<NB; ++k) for(i=0; i<TAILLE; ++i)
        GFq.axpyin(z3, z1, z2);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " MulAdd IN place: " << 2*coef / tim.usertime()
         << "Mop/s\n" << tim << ", ex: " << z3 << std::endl << std::flush;

    tim.clear();tim.start();
    for (k=0; k<NB; ++k) for(i=0; i<TAILLE; ++i)
        GFq.axpy(z3, z1, z2, z23);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " MulAdd: " << 2*coef / tim.usertime()
         << "Mop/s\n" << tim << ", ex: " << z3 << std::endl << std::flush;

    tim.clear();tim.start();
    for (k=0; k<NB; ++k) for(i=0; i<TAILLE; ++i)
        GFq.neg(z3, z1);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " Neg: " << coef / tim.usertime()
         << "Mop/s\n" << tim << ", ex: " << z3 << std::endl << std::flush;

    return 0;
}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

