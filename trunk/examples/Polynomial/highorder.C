// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software. 
// see the COPYRIGHT file for more details.

#include <iostream>
#include <givaro/givrandom.h>
#include <givaro/givtimer.h>
#include <givaro/givgfq.h>
#include <givaro/givpoly1.h>
#include "givtruncdomain.h"
#include "givhighorder.h"


typedef GFqDom<int> Field;
typedef Poly1Dom< Field, Dense> Polys;
typedef HighOrder< Field > HighOrders;


long long count = 0;
Timer Ttaylor, Thighorder; 



bool TestFracDevel(const HighOrders& HO101, const Polys::Element P, const Polys::Element oQ, Degree a, Degree b) {
    ++count;
//     std::cerr << "------------------------------------------------------------" << count <<std::endl;

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
    
    bool success = HO101.gettruncdom().areEqual(TTy17,TTy17b);

    if (! success) {
        std::cerr << "---------- BEG ERROR ----------" <<std::endl;
        HO101.getpoldom().write(std::cerr << "P:=", P) << ';' << std::endl;
        HO101.getpoldom().write(std::cerr << "Q:=", Q) << ';' << std::endl;
        HO101.getpoldom().write(std::cerr << "TTy:=", TTy) << ';' << std::endl;

        Fra._num = HO101.getpoldom().one;
        Polys::Element Tay;
        HO101.taylor(Tay, Fra, b);
        HO101.getpoldom().write(std::cerr << "Tay:=", Tay) << ';' << std::endl;
        

        HO101.write(std::cerr << "TTy17:=", TTy17) << ';' << std::endl;
        HO101.write(std::cerr << "TTy17b:=", TTy17b) << ';' << std::endl;
        std::cerr << "---------- END ERROR ----------" <<std::endl;
    } 
        

    return success;
}
    





int main(int argc, char ** argv) {
    long seed = (argc>1?atoi(argv[1]):BaseTimer::seed());
    std::cerr << "seed: " << seed << std::endl;

    Field Z101( 101, 1 );  // integers modulo 101

    // Polynomials over Z101, with X as indeterminate
    Polys DP101( Z101, Indeter("X") );
    Polys::Element P, Q, R, monomial;

    GivRandom generator(seed);
    long deg1 = generator() % 6;
    long deg2 = generator() % 6;
    long deg3 = generator() % 155;
    long v1 = generator() % 195;
    long v2 = v1 + (generator() % 5);
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
    
    size_t e;

    HighOrders::Truncated G0;
    HO101.GammaId(G0, S, dS, e, P, dP);
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

    success &= TestFracDevel(HO101, P, Q, 100-17,100);
    success &= TestFracDevel(HO101, P, Q, 100-27,100);
    success &= TestFracDevel(HO101, P, Q, 100-87,100);
    success &= TestFracDevel(HO101, P, Q, 0,100);
    success &= TestFracDevel(HO101, P, Q, 1,100);
    success &= TestFracDevel(HO101, P, Q, 2,100);
    success &= TestFracDevel(HO101, P, Q, 100,100);
    success &= TestFracDevel(HO101, P, Q, 99,100);
    success &= TestFracDevel(HO101, P, Q, 98,100);


// Error 857475
//     for(size_t i=0; i<452; ++i) {
//         long deg1 = generator() % 16;
//         long deg2 = generator() % 15;
//         long v1 = generator() % 29;
//         long v2 = v1 + (generator() % 24);
//         DP101.random(generator, P, Degree(deg1) );
//         DP101.random(generator, Q, Degree(deg2) );
//         success &= TestFracDevel(HO101, P, Q, v1, v2);
//     }

// Error 201149
    for(size_t i=0; i<2000; ++i) {
        long deg1 = generator() % 66;
        long deg2 = generator() % 65;
        long v1 = generator() % 19195;
        long v2 = v1 + (generator() % 45);
        DP101.random(generator, P, Degree(deg1) );
        DP101.random(generator, Q, Degree(deg2) );
        success &= TestFracDevel(HO101, P, Q, v1, v2);
    }

    std::cerr << "Taylor: " << Ttaylor << std::endl;
    std::cerr << "HighOrder: " << Thighorder << std::endl;
    


    if (! success) {
        std::cerr << "Error: " << seed << std::endl;
    } else {
        std::cerr << "Success:" << count << std::endl;
    }

    return 0;
}
