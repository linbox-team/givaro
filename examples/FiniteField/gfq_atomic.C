#include <stdio.h>
#include <stdlib.h>
#include <givaro/givtimer.h>
#include <givaro/givrandom.h>
#include <givaro/givgfq.h>

#ifndef GIVMIN
#define GIVMIN(a,b) (((a)<(b))?(a):(b))
#endif

typedef GFqDom<long> Domain;
typedef GFqDom<long>::element Modulo;

// NB: number of iterations, TAILLE: vector size
#ifndef NB
#define NB 500000
#endif
#ifndef TAILLE
#define TAILLE 256
#endif

int main(int argc, char ** argv) {

    GivRandom generator;
    long P (65521);     // argv[1] : characteristic
    long expo(1);       // argv[2] : exponent
    int offset = 0;
    if (argc > ++offset) P = atoi( argv[ offset ] );
    if (argc > ++offset) expo = atoi( argv[ offset ] );
				        
    Timer inver;
    inver.clear();
    inver.start();
    Domain GFq(P, expo);  // Buiding of finite field with P^expo elements
							    
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
         << "Mop/s\n" << tim << std::endl << std::flush;

    tim.clear();tim.start();
    for (k=0; k<NB; ++k)  for(i=0; i<TAILLE; ++i)
        GFq.add(z3, z1, z2);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " Add: " << coef / tim.usertime() 
         << "Mop/s\n" << tim << std::endl << std::flush;

    tim.clear();tim.start();
    for (k=0; k<NB; ++k)  for(i=0; i<TAILLE; ++i)
        GFq.axpyin(z3, z1, z2);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " MulAdd IN place: " << 2*coef / tim.usertime() 
         << "Mop/s\n" << tim << std::endl << std::flush;

    tim.clear();tim.start();
    for (k=0; k<NB; ++k)  for(i=0; i<TAILLE; ++i)
        GFq.axpy(z3, z1, z2, z23);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " MulAdd: " << 2*coef / tim.usertime() 
         << "Mop/s\n" << tim << std::endl << std::flush;

    tim.clear();tim.start();
    for (k=0; k<NB; ++k)  for(i=0; i<TAILLE; ++i)
        GFq.sub(z3, z1, z2);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " Sub: " << coef / tim.usertime() 
         << "Mop/s\n" << tim << std::endl << std::flush;
/*
    tim.clear();tim.start();
    for (k=0; k<NB; ++k)  for(i=0; i<TAILLE; ++i)
        GFq.div(z3, z1, z2);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " Div: " << coef / tim.usertime() 
         << "Mop/s\n" << tim << std::endl << std::flush;
*/

    tim.clear();tim.start();
    for (k=0; k<NB; ++k) for(i=0; i<TAILLE; ++i)
        GFq.neg(z3, z1);
    std::cout << NB << " * " << TAILLE << " Neg: " << coef / tim.usertime() 
         << "Mop/s\n" << tim << std::endl << std::flush;

            GFq.random(generator,z1) ;
            GFq.nonzerorandom(generator,z2) ;
            GFq.random(generator,z23) ;


    tim.clear();tim.start();
    for (k=0; k<NB; ++k) for(i=0; i<TAILLE; ++i)
        GFq.add(z3, z1, z2);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " Add: " << coef / tim.usertime() 
         << "Mop/s\n" << tim << std::endl << std::flush;

    tim.clear();tim.start();
    for (k=0; k<NB; ++k) for(i=0; i<TAILLE; ++i)
        GFq.sub(z3, z1, z2);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " Sub: " << coef / tim.usertime() 
         << "Mop/s\n" << tim << std::endl << std::flush;

    tim.clear();tim.start();
    for (k=0; k<NB; ++k) for(i=0; i<TAILLE; ++i)
        GFq.axpyin(z3, z1, z2);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " MulAdd IN place: " << 2*coef / tim.usertime()
         << "Mop/s\n" << tim << std::endl << std::flush;

    tim.clear();tim.start();
    for (k=0; k<NB; ++k) for(i=0; i<TAILLE; ++i)
        GFq.axpy(z3, z1, z2, z23);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " MulAdd: " << 2*coef / tim.usertime() 
         << "Mop/s\n" << tim << std::endl << std::flush;

    tim.clear();tim.start();
    for (k=0; k<NB; ++k) for(i=0; i<TAILLE; ++i)
        GFq.neg(z3, z1);
    tim.stop();
    std::cout << NB << " * " << TAILLE << " Neg: " << coef / tim.usertime() 
         << "Mop/s\n" << tim << std::endl << std::flush;

    return 0;
}

