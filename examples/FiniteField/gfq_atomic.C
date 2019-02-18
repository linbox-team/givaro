// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

/*! @file examples/FiniteField/gfq_atomic.C
 * @ingroup examples
 * @ingroup finitefields
 * @example examples/FiniteField/gfq_atomic.C
 * @brief NO DOC
 */
#include <stdio.h>
#include <stdlib.h>
#include <givaro/givtimer.h>
#include <givaro/givrandom.h>
#include <givaro/gfq.h>

using namespace Givaro;



#ifndef GIVMIN
#define GIVMIN(a,b) (((a)<(b))?(a):(b))
#endif

typedef GFqDom<int64_t> Domain;
typedef GFqDom<int64_t>::Element Modulo;

// NB: number of iterations, TAILLE: vector size
#ifndef NB
#define NB 1000000
#endif
#ifndef TAILLE
#define TAILLE 256
#endif

int main(int argc, char ** argv)
{

    GivRandom generator;
    int64_t P (65521);     // argv[1] : characteristic
    int64_t expo(1);       // argv[2] : exponent
    int offset = 0;
    if (argc > ++offset) P = atoi( argv[ offset ] );
    if (argc > ++offset) expo = atoi( argv[ offset ] );

    Timer inver;
    inver.clear();
    inver.start();
    Domain GFq((uint64_t)P, (uint64_t)expo);  // Buiding of finite field with P^expo Elements

    Modulo * z1 = new Modulo[TAILLE], * z2 = new Modulo[TAILLE], * z23 = new Modulo[TAILLE], * z3 = new Modulo[TAILLE];

    //    int64_t seuil = GIVMIN(P*2,TAILLE);

    std::cout << "."<< std::flush;

    for(int i=0; i<TAILLE; ++i){
        GFq.random(generator,z1[i]) ;
        GFq.nonzerorandom(generator,z2[i]) ;
        GFq.random(generator,z23[i]) ;
    }


    std::cout << "."<< std::flush;
    inver.stop();

    std::cout << "Init et Tableaux des  modulo " << P << " :\n" << inver << std::endl << std::flush;

    double coef = (double)NB*TAILLE / 1e6;

    Timer tim;tim.clear();tim.start();
    for (int k=0; k<NB; ++k)  for(int i=0; i<TAILLE; ++i)
        GFq.mul(z3[i], z1[i], z2[i]);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " Mul: " << coef / tim.usertime()
    << "Mop/s\n" << tim << ", ex: " << z3[0] << std::endl << std::flush;

    tim.clear();tim.start();
    for (int k=0; k<NB; ++k)  for(int i=0; i<TAILLE; ++i)
        GFq.add(z3[i], z1[i], z2[i]);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " Add: " << coef / tim.usertime()
    << "Mop/s\n" << tim << ", ex: " << z3[0] << std::endl << std::flush;

    tim.clear();tim.start();
    for (int k=0; k<NB; ++k)  for(int i=0; i<TAILLE; ++i)
        GFq.axpyin(z3[i], z1[i], z2[i]);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " MulAdd IN place: " << 2*coef / tim.usertime()
    << "Mop/s\n" << tim << ", ex: " << z3[0] << std::endl << std::flush;

    tim.clear();tim.start();
    for (int k=0; k<NB; ++k)  for(int i=0; i<TAILLE; ++i)
        GFq.axpy(z3[i], z1[i], z2[i], z23[i]);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " MulAdd: " << 2*coef / tim.usertime()
    << "Mop/s\n" << tim << ", ex: " << z3[0] << std::endl << std::flush;

    tim.clear();tim.start();
    for (int k=0; k<NB; ++k)  for(int i=0; i<TAILLE; ++i)
        GFq.sub(z3[i], z1[i], z2[i]);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " Sub: " << coef / tim.usertime()
    << "Mop/s\n" << tim << ", ex: " << z3[0] << std::endl << std::flush;
#if 0
    tim.clear();tim.start();
    for (int k=0; k<NB; ++k)  for(int i=0; i<TAILLE; ++i)
        GFq.div(z3[i], z1[i], z2[i]);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " Div: " << coef / tim.usertime()
    << "Mop/s\n" << tim << ", ex: " << z3[0] << std::endl << std::flush;
#endif

    tim.clear();tim.start();
    for (int k=0; k<NB; ++k) for(int i=0; i<TAILLE; ++i)
        GFq.neg(z3[i], z1[i]);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " Neg: " << (coef / tim.usertime())
    << "Mop/s\n" << tim << ", ex: " << z3[0] << std::endl << std::flush;

    for(int i=0; i<TAILLE; ++i){
        GFq.random(generator,z1[i]) ;
        GFq.nonzerorandom(generator,z2[i]) ;
        GFq.random(generator,z23[i]) ;
    }


    tim.clear();tim.start();
    for (int k=0; k<NB; ++k) for(int i=0; i<TAILLE; ++i)
        GFq.add(z3[i], z1[i], z2[i]);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " Add: " << coef / tim.usertime()
    << "Mop/s\n" << tim << ", ex: " << z3[0] << std::endl << std::flush;

    tim.clear();tim.start();
    for (int k=0; k<NB; ++k) for(int i=0; i<TAILLE; ++i)
        GFq.sub(z3[i], z1[i], z2[i]);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " Sub: " << coef / tim.usertime()
    << "Mop/s\n" << tim << ", ex: " << z3[0] << std::endl << std::flush;

    tim.clear();tim.start();
    for (int k=0; k<NB; ++k) for(int i=0; i<TAILLE; ++i)
        GFq.axpyin(z3[i], z1[i], z2[i]);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " MulAdd IN place: " << 2*coef / tim.usertime()
    << "Mop/s\n" << tim << ", ex: " << z3[0] << std::endl << std::flush;

    tim.clear();tim.start();
    for (int k=0; k<NB; ++k) for(int i=0; i<TAILLE; ++i)
        GFq.axpy(z3[i], z1[i], z2[i], z23[i]);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " MulAdd: " << 2*coef / tim.usertime()
    << "Mop/s\n" << tim << ", ex: " << z3[0] << std::endl << std::flush;

    tim.clear();tim.start();
    for (int k=0; k<NB; ++k) for(int i=0; i<TAILLE; ++i)
        GFq.neg(z3[i], z1[i]);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " Neg: " << coef / tim.usertime()
    << "Mop/s\n" << tim << ", ex: " << z3[0] << std::endl << std::flush;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
