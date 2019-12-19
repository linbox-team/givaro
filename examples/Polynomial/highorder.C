// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
/*! @file examples/Polynomial/highorder.C
 * @ingroup examples
 * @ingroup polynomials
 * @example examples/Polynomial/highorder.C
 * @brief NO DOC
 */

#include <iostream>
#ifndef __GIVARO_DEBUG
#define __GIVARO_DEBUG 1
#endif
#include <givaro/givrandom.h>
#include <givaro/givtimer.h>
#include <givaro/gfq.h>
#include <givaro/givpoly1.h>
#include <givaro/givtruncdomain.h>
#include <givaro/givhighorder.h>

using namespace Givaro;




typedef GFqDom<int> Field;
typedef Poly1Dom< Field, Dense> Polys;
typedef HighOrder< Field > HighOrders;


long long TFDcount = 0;
long long FiDcount = 0;
Timer Ttaylor, Thighorder, Tfiduccia;



bool TestFracDevel(const HighOrders& HO101, const Polys::Element P, const Polys::Element oQ, Degree a, Degree b) {
    bool success = false, successF=false;

    try {


        ++TFDcount;
        //     std::cerr << "------------------------------------------------------------" << TFDcount <<std::endl;

        HighOrders::Truncated TQ;
        Degree dp; HO101.getpoldom().degree(dp, P);
        HO101.gettruncdom().assign(TQ, oQ,0,dp-1);
        Polys::Element Q;
        HO101.gettruncdom().convert(Q, TQ);
        //HO101.write(std::cout << "Q:=", TQ) << ';' << std::endl;

        HighOrders::Element Fra;
        Fra._num = Q;
        Fra._den = P;

        Polys::Element TTy;
        Timer tt; tt.clear(); tt.start();
        HO101.taylor(TTy, Fra, b);
        tt.stop();
        Ttaylor+=tt;

        HighOrders::Truncated TTy17;
        HO101.gettruncdom().assign(TTy17, TTy, a, b);

        //HO101.write(std::cout << "TTy17:=", TTy17) << ';' << std::endl;


        HighOrders::Truncated TTy17b;

        tt.clear(); tt.start();
        HO101.FracDevel(TTy17b, Fra, a, b);
        tt.stop();
        Thighorder+=tt;


        //HO101.write(std::cout << "TTy17b:=", TTy17b) << ';' << std::endl;

        success = HO101.gettruncdom().areEqual(TTy17,TTy17b);

        if (! success) {
            std::cerr << "---------- BEG ERROR ----------" <<std::endl;
            HO101.getpoldom().write(std::cerr << "Q:=", Q) << ';' << std::endl;
            HO101.getpoldom().write(std::cerr << "P:=", P) << ';' << std::endl;
            HO101.getpoldom().write(std::cerr << "TTy:=", TTy) << ';' << std::endl;

            //         Fra._num = HO101.getpoldom().one;
            //         Polys::Element Tay;
            //         HO101.taylor(Tay, Fra, b);
            //         HO101.getpoldom().write(std::cerr << "Tay:=", Tay) << ';' << std::endl;


            HO101.write(std::cerr << "TTy17:=", TTy17) << ';' << std::endl;
            HO101.write(std::cerr << "TTy17b:=", TTy17b) << ';' << std::endl;
            std::cerr << "---------- END ERROR ----------" <<std::endl;
        }


        HighOrders::Truncated TTy17c; HO101.gettruncdom().init(TTy17c);
        HighOrders::Truncated TTy17d; HO101.gettruncdom().init(TTy17d);

        tt.clear(); tt.start();
        //     HO101.Fiduccia(TTy17c, Fra, b);
        //     for(Degree d=a; d<b; ++d) {
        //         HO101.Fiduccia(TTy17d, Fra, d);
        //         HO101.gettruncdom().addin(TTy17c, TTy17d);
        //     }
        HO101.Fiduccia(TTy17c, Fra, a, b);
        tt.stop();
        Tfiduccia+=tt;


        //HO101.write(std::cout << "TTy17b:=", TTy17b) << ';' << std::endl;

        successF = HO101.gettruncdom().areEqual(TTy17,TTy17c);

        if (! successF) {
            std::cerr << "---------- BEG ERROR ----------" <<std::endl;
            HO101.getpoldom().write(std::cerr << "Q:=", Q) << ';' << std::endl;
            HO101.getpoldom().write(std::cerr << "P:=", P) << ';' << std::endl;
            HO101.getpoldom().write(std::cerr << "TTy:=", TTy) << ';' << std::endl;

            //         Fra._num = HO101.getpoldom().one;
            //         Polys::Element Tay;
            //         HO101.taylor(Tay, Fra, b);
            //         HO101.getpoldom().write(std::cerr << "Tay:=", Tay) << ';' << std::endl;


            HO101.write(std::cerr << "TTy17:=", TTy17) << ';' << std::endl;
            HO101.write(std::cerr << "TTy17c:=", TTy17c) << ';' << std::endl;
            std::cerr << "---------- END ERROR ----------" <<std::endl;
        } else {
            ++FiDcount;
        }


    } catch (GivError &e) {
        std::cerr << e << std::endl;
    }

    return success && successF;
}






int main(int argc, char ** argv) {


    long numb = (argc>1?atoi(argv[1]):200);
    std::cerr << "numb: " << numb << std::endl;
    long tttn = (argc>2?atoi(argv[2]):100);
    std::cerr << "tttn: " << tttn << std::endl;
    long seed = (argc>3?atoi(argv[3]):BaseTimer::seed());
    std::cerr << "seed: " << seed << std::endl;

    Field Z101( 101, 1 );  // integers modulo 101

    // Polynomials over Z101, with X as indeterminate
    Polys DP101( Z101, Indeter("X") );
    Polys::Element P, Q, R, monomial;

    GivRandom generator((unsigned long)seed);
    long deg1 = (long) generator() % 6;
    long deg2 = (long) generator() % 6;
    long deg3 = (long) generator() % 155;
    // long v1 = generator() % 195;
    // long v2 = v1 +  (generator() % 5);
    DP101.random(generator, P, Degree(deg1) );
    DP101.random(generator, Q, Degree(deg2) );
    DP101.random(generator, R, Degree(deg3) );
    Degree dP; DP101.degree(dP,P);
    Degree dQ; DP101.degree(dQ,Q);
    Degree dR; DP101.degree(dR,R);

    DP101.write(std::cout << "P:=", P) << ';' << std::endl;

    HighOrders HO101( Z101, Indeter("X") );

    HighOrders::Element F;
    F._num = DP101.one;
    F._den = P;

    Polys::Element Tay;

    HO101.taylor(Tay, F, 128);
    DP101.write(std::cout << "Tay:=", Tay) << ';' << std::endl;


    Polys::Element S; Degree dS;

    size_t e = 0 ; // initialisé à quoi ? dans GammaId, k0 est const...

    HighOrders::Truncated G0;
    HO101.GammaId(G0, S, dS, (long)e, P, dP);
    std::cout << "e:=" << e << ';' << std::endl;
    HO101.write(std::cout << "G0:=", G0) << ';' << std::endl;
    DP101.write(std::cout << "S:=", S) << ';' << std::endl;

    std::vector<HighOrders::Truncated> Gam, T;
    std::vector<Degree> D;

    Polys::Element nTay; Degree dT;

    HO101.highorder(Gam, T, D, nTay, dT, Degree(970), Degree(1000), P, dP);
    HO101.write(std::cout << "Gam:=", Gam.back()) << ';' << std::endl;

    HighOrders::Truncated TP; HO101.gettruncdom().assign(TP, P);

    HighOrders::Truncated I;
    HO101.Inverse(I, 113, 127, TP, dP, nTay, dT, Gam, T, D);

    HO101.write(std::cout << "I]_155^175:=", I) << ';' << std::endl;

    Ttaylor.clear();
    Thighorder.clear();

    bool success=true;

    success &= TestFracDevel(HO101, P, Q, tttn-(17*tttn)/100,tttn);
    success &= TestFracDevel(HO101, P, Q, tttn-(27*tttn)/100,tttn);
    success &= TestFracDevel(HO101, P, Q, tttn-(87*tttn)/100,tttn);
    success &= TestFracDevel(HO101, P, Q, 0,tttn);
    success &= TestFracDevel(HO101, P, Q, 1,tttn);
    success &= TestFracDevel(HO101, P, Q, 2,tttn);
    success &= TestFracDevel(HO101, P, Q, tttn,tttn);
    success &= TestFracDevel(HO101, P, Q, tttn-1,tttn);
    success &= TestFracDevel(HO101, P, Q, tttn-2,tttn);


    for(long i=0; i<numb; ++i) {
        long Deg1 = (long)generator() % ((66*tttn)/100);
        long Deg2 = (long)generator() % ((65*tttn)/100);
        long v1 = (long)generator() % ((19195*tttn)/100);
        long v2 = v1 + ((long)generator() % (long)((45*tttn)/100));
        DP101.random(generator, P, Degree(Deg1) );
        DP101.random(generator, Q, Degree(Deg2) );
        success &= TestFracDevel(HO101, P, Q, v1, v2);
    }

    std::cerr << "Taylor: " << Ttaylor << std::endl;
    std::cerr << "HighOrder: " << Thighorder << std::endl;
    std::cerr << "Fiduccia: " << Tfiduccia << std::endl;



    if (! success) {
        std::cerr << "Error: " << seed << std::endl;
    } else {
        std::cerr << "Success:" << TFDcount << std::endl;
    }
    std::cerr << "Success F:" << FiDcount << std::endl;

    return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
