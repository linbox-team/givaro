// =================================================================== //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Contributor: Jack Dubrois < Jacques.Dubrois@imag.fr>
// Time-stamp: <12 Jun 15 18:24:19 Jean-Guillaume.Dumas@imag.fr>
//
// Primality check using Probabilistic Lucas    /////////////////////////
// i.e. Primitive Root with choosen probability /////////////////////////
// =================================================================== //

/*! @file examples/Integer/ProbLucas.C
 * @ingroup examples
 * @ingroup integers
 * @example examples/Integer/ProbLucas.C
 * @brief NO DOC
 */
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
using namespace std;
#include <stdlib.h>
#include <cmath>
#include <givaro/givintprime.h>
#include <givaro/givintnumtheo.h>
#include <givaro/givtimer.h>
#include <givaro/givinit.h>
#include <givaro/modular-integer.h>
#include <givaro/givinteger.h>
#include <givaro/givrandom.h>
#include <givaro/givinteger.h>


using namespace Givaro;





IntFactorDom<> IP;
#define GIVABSDIV(a,b) ((a)<(b)?((a)/(b)):((b)/(a)))
#ifndef GIVMIN
#define GIVMIN(a,b) ((a)<(b)?(a):(b))
#endif
#define ProbLucas_factor_first_primes(tmp,n) (tmp = IP.isZero(IP.mod(tmp,n,23))?23:( IP.isZero(IP.mod(tmp,n,19))?19:( IP.isZero(IP.mod(tmp,n,17))?17:  (IP.isZero(IP.mod(tmp,n,2))?2:( IP.isZero(IP.mod(tmp,n,3))?3:( IP.isZero(IP.mod(tmp,n,5))?5:( IP.isZero(IP.mod(tmp,n,7))?7: ( IP.isZero(IP.mod(tmp,n,11))?11:13 ))))))))

#define ProbLucas_factor_second_primes(tmp,n) (tmp = IP.isZero(IP.mod(tmp,n,31))?31:( IP.isZero(IP.mod(tmp,n,29))?29: ( IP.isZero(IP.mod(tmp,n,37))?37: ( IP.isZero(IP.mod(tmp,n,41))?41:( IP.isZero(IP.mod(tmp,n,43))?43:  ( IP.isZero(IP.mod(tmp,n,71))?71:( IP.isZero(IP.mod(tmp,n,67))?67:( IP.isZero(IP.mod(tmp,n,61))?61:( IP.isZero(IP.mod(tmp,n,59))?59: ( IP.isZero(IP.mod(tmp,n,53))?53:( IP.isZero(IP.mod(tmp,n,47))?47: ( IP.isZero(IP.mod(tmp,n,97))?97: ( IP.isZero(IP.mod(tmp,n,89))?89:( IP.isZero(IP.mod(tmp,n,83))?83:( IP.isZero(IP.mod(tmp,n,79))?79:73)))))))))))))))



Integer& MyPollard(GivRandom& gen, Integer& g, const Integer& n, const unsigned long threshold)
{
    g=1;
    Integer m(0U), x, y, t;
    static Integer PROD_first_primes(223092870);
    static Integer PROD_second_primes("10334565887047481278774629361");
    if (isOne(gcd(y,n,PROD_first_primes))) {
        if (isOne(gcd(y,n,PROD_second_primes))) {

            IP.random(gen, y, n);
            unsigned long p(1);
            for(unsigned long c = 0; isOne(g) && (++c < threshold); ) {
                if(  p == c ) {
                    x=y;
                    p <<= 1;
                }
                // Pollard fctin
                IP.mulin(y,y);
                IP.addin(y,1U);
                IP.modin(y,n);
                gcd(g,IP.sub(t,y,x),n);
            }
            return g;
        } else {
            return ProbLucas_factor_second_primes(g,n);
        }
    } else {
        return ProbLucas_factor_first_primes(g,n);
    }

}







unsigned long Revert(const Integer p, const double epsilon, double firstguess)
{
    // unsigned long L;
    double t1, t4, t8, t34(firstguess), dL;
    t1 = (double)p-1.0;
    t4 = 2.0/t1+1.0;
    t8 = logtwo(p)*log(2.0)-log(2.0);            // log( p/2 )
    do {
        double t3, t5, t7, t9, t11, t12, t18, t22, t23 ;
        //          std::cerr << "dL: " << t34 << std::endl;
        dL = t34;
        t3 = 1.0/dL;
        t5 = t3*t3;
        t7 = 1.0-t5;
        t9 = 2.0*log(dL);
        t11 = t8/t9;
        t12 = ::pow(t7,t11);
        t18 = t9*t9;
        t22 = log(t7);
        t23 = (-t8/t18*t3*t22+t11*t5*t3/t7);
        t34 = dL-(  1.0 - (1.0-epsilon)/(t4*t12))/(2.0*t23);
    } while( (GIVABSDIV(t34,dL) < 0.95) && (dL<1048576.0)  );
    return /* L =*/ (unsigned long)GIVMIN(dL,1048576.0);
}

unsigned long Revert(const Integer p, const double epsilon)
{
    return Revert(p, epsilon, ::sqrt(1.2/epsilon));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Recherche de la racine primitive ///////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool ProbLucas(const Integer n, const double orig_epsilon)
{
#ifdef __GMP_PLUSPLUS__
    Integer::seeding( (unsigned long)BaseTimer::seed() );
#endif
    GivRandom generator;

    Integer Q=n-1,a,q,tmp(1);
    Integer nmu=Q;

    double P = 1.0;
    double epsilon = orig_epsilon;

    // Miller-Rabin
    unsigned long s=0;
    for( ; !( (int)Q & 0x1) ; Q>>=1, ++s) { }
    for(unsigned int i=0; ; ) {
        IP.random(generator,a,n);
        IP.powmod(tmp,a,Q,n);
        if (tmp!=1) {
            if (tmp != nmu) {
                for(i=(unsigned int) s; --i>0;) {
                    tmp *= tmp;
                    tmp %= n;
                    if ((tmp == 1) || (tmp == nmu))
                        break;
                }
                if (tmp != nmu)
                    return 0;
                if (i == 0)
                    break;
            }
            else {
                if (s == 1) break;
            }
        }
        epsilon *= 4.0;
        P /= 2.0;
    }

    unsigned long L = Revert(Q,epsilon);
    Integer trn, r, b, expo; root( trn, n, 3); trn *= trn;
    q = 2;
    expo = nmu;
    MyPollard(generator, q, Q, L);
    if (q == 1) {
        q = Q;
        expo = nmu; expo /= q;
        IP.random(generator, a, n);
        IP.powmod(tmp,a,expo,n);
        if (tmp == 1) {
            P /= (double)(q);
            if (P < epsilon) {
                std::cerr << "Composite with probability > 1-" << P << std::endl;
                return 0;
            }
        }
        Q = 1;
    } else {
        expo = nmu; expo /= q;
        r=0;
        Integer::divexact(b, Q, q);
        while( isZero(r) ) {
            Q.copy(b);
            Integer::divmod( b, r, Q, q );
        }
        L = Revert(Q,epsilon, (double)L);
    }

    while ( Q > trn ) {

        IP.random(generator, a, n);
        IP.powmod(tmp,a,expo,n);
        if (tmp == 1) {
            P /= (double)(q);
            if (P < epsilon) {
                std::cerr << "Composite with probability > 1-" << P << std::endl;
                return 0;
            }
        } else {
            if (IP.gcd(q,(tmp-1U),n) != 1) {
                std::cerr << "Factor found : " << q << std::endl;
                return 0;
            }
            MyPollard(generator, q, Q, L);
            if (q == 1) {
                q = Q;
                expo = nmu; expo /= q;
                IP.random(generator, a, n);
                IP.powmod(tmp,a,expo,n);
                if (tmp == 1) {
                    P /= (double)(q);
                    if (P < epsilon) {
                        std::cerr << "Composite with probability > 1-" << P << std::endl;
                        return 0;
                    }
                }
                break;
            } else {
                expo = nmu; expo /= q;
                r=0;
                Integer::divexact(b, Q, q);
                while( isZero(r) ) {
                    Q.copy(b);
                    Integer::divmod( b, r, Q, q );
                }
                L = Revert(Q,epsilon, (double)L);
            }
        }
    }
    if (! IP.isprime(q))
        std::cerr << "Prime with probability > 1-" << orig_epsilon << std::endl;
    return 1;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int main (int argc, char * * argv)
{

    Integer P;
    if (argc > 1) P = Integer(argv[1]); else std::cin >> P;
    double epsilon = argc > 2 ? atof(argv[2]) : 0.000001;
    unsigned int NB = argc > 3 ? (unsigned int)atoi(argv[3]) : 1U;

    //    std::cerr << "P: " << P << " ; proba: " << epsilon << std::endl;

    bool a1(true);
    Timer tim; tim.clear();
    for(unsigned int i = 0; i < NB; ++i) {
        Timer tt; tt.clear(); tt.start();
        a1 = ProbLucas(P, epsilon);
        tt.stop();
        //        std::cout << a1 << std::endl;
        //     	  std::cerr << tt << std::endl;
        tim += tt;
    }
    std::cout << (a1?"prime":"composite") << endl;
    std::cerr << tim << std::endl;
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
