// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

#include <iostream>
#include <givaro/gfq.h>
#include <givaro/givrandom.h>
#include <givaro/givtimer.h>
#include <givaro/givpoly1.h>
#include <givaro/givtruncdomain.h>

using namespace Givaro;

long long TTcount = 0;

bool TestAdd(const TruncDom< GFqDom<int> >& DP, const TruncDom< GFqDom<int> >::Element& P, const TruncDom< GFqDom<int> >::Element& Q, size_t d1, size_t d2)
{
    ++TTcount;
    TruncDom< GFqDom<int> >::Element R, T, V, S, U, W;
    DP.add ( R, P, Q, (long)d1, (long)d2); // R = P+Q;
    //     DP.write( DP.write(
    //         std::cout << "[(" , P ) << ") + (", Q) << ")]_" << d1 << '^' << d2 ;
    //     DP.write(std::cout << " = " , R) << std::endl;

    DP.add (T, P, Q);
    V=T;
    DP.truncin(V,(long)d1,(long)d2);

    S = P; U = P;
    DP.addin(S,Q,(long)d1,(long)d2);
    DP.addin(U,Q);
    W=U;
    DP.truncin(W,(long)d1,(long)d2);

    if( DP.areNEqual( V, R) || DP.areNEqual(W, S) || DP.areNEqual(V, W) ) {
        std::cerr << "ERROR ADD:" << TTcount << std::endl;
        DP.write(std::cout << "  R: " , R) << std::endl;
        DP.write(std::cout << "  T: " , T) << std::endl;
        DP.write(std::cout << "  V: " , V) << std::endl;
        Degree dP; DP.degree(dP,P);
        Degree dQ; DP.degree(dQ,Q);
        Degree vP; DP.val(vP,P);
        Degree vQ; DP.val(vQ,Q);
        std::cerr << "vR: " << vP << ", dR: " << dP << ", vP: " << vQ << ", dP: " << dQ << ", v: " << d1 << ", d: " << d2 << std::endl;
        return false;
    }
    return true;
}

bool TestSub(const TruncDom< GFqDom<int> >& DP, const TruncDom< GFqDom<int> >::Element& P, const TruncDom< GFqDom<int> >::Element& Q, size_t d1, size_t d2)
{
    ++TTcount;
    TruncDom< GFqDom<int> >::Element R, T, V, S, U, W;
    DP.sub ( R, P, Q, (long)d1, (long)d2); // R = P-Q;
    //     DP.write( DP.write(
    //         std::cout << "[(" , P ) << ") + (", Q) << ")]_" << d1 << '^' << d2 ;
    //     DP.write(std::cout << " = " , R) << std::endl;

    DP.sub (T, P, Q);
    V=T;
    DP.truncin(V,(long)d1,(long)d2);

    S = P; U = P;
    DP.subin(S,Q,(long)d1,(long)d2);
    DP.subin(U,Q);
    W=U;
    DP.truncin(W,(long)d1,(long)d2);

    if( DP.areNEqual( V, R) || DP.areNEqual(W, S) || DP.areNEqual(V, W) ) {
        std::cerr << "ERROR SUB:" << TTcount << std::endl;
        DP.write(std::cout << "  R: " , R) << std::endl;
        DP.write(std::cout << "  T: " , T) << std::endl;
        DP.write(std::cout << "  V: " , V) << std::endl;
        Degree dP; DP.degree(dP,P);
        Degree dQ; DP.degree(dQ,Q);
        Degree vP; DP.val(vP,P);
        Degree vQ; DP.val(vQ,Q);
        std::cerr << "vR: " << vP << ", dR: " << dP << ", vP: " << vQ << ", dP: " << dQ << ", v: " << d1 << ", d: " << d2 << std::endl;
        return false;
    }
    return true;
}

bool TestMul(const TruncDom< GFqDom<int> >& DP, const TruncDom< GFqDom<int> >::Element& P, const TruncDom< GFqDom<int> >::Element& Q, size_t d1, size_t d2)
{
    ++TTcount;
    TruncDom< GFqDom<int> >::Element R, T, V, S, U, W;
    DP.mul ( R, P, Q, (long)d1, (long)d2); // R = P*Q;
    //     DP.write( DP.write(
    //         std::cout << "[(" , P ) << ") + (", Q) << ")]_" << d1 << '^' << d2 ;
    //     DP.write(std::cout << " = " , R) << std::endl;

    DP.mul (T, P, Q);
    V=T;
    DP.truncin(V,(long)d1,(long)d2);

    S = P; U = P;
    DP.mulin(S,Q,(long)d1,(long)d2);
    DP.mulin(U,Q);
    W=U;
    DP.truncin(W,(long)d1,(long)d2);

    if( DP.areNEqual( V, R) || DP.areNEqual(W, S) || DP.areNEqual(V, W) ) {
        std::cerr << "ERROR MUL:" << TTcount << std::endl;
        DP.write(std::cout << "  R: " , R) << std::endl;
        DP.write(std::cout << "  T: " , T) << std::endl;
        DP.write(std::cout << "  V: " , V) << std::endl;
        Degree dP; DP.degree(dP,P);
        Degree dQ; DP.degree(dQ,Q);
        Degree vP; DP.val(vP,P);
        Degree vQ; DP.val(vQ,Q);
        std::cerr << "vR: " << vP << ", dR: " << dP << ", vP: " << vQ << ", dP: " << dQ << ", v: " << d1 << ", d: " << d2 << std::endl;
        return false;
    }
    return true;
}

bool TestAxpy(const TruncDom< GFqDom<int> >& DP, const TruncDom< GFqDom<int> >::Element& P, const TruncDom< GFqDom<int> >::Element& Q, const TruncDom< GFqDom<int> >::Element& G, size_t d1, size_t d2)
{
    ++TTcount;
    TruncDom< GFqDom<int> >::Element R, T, V, S, U, W;
    DP.axpy ( R, P, Q, G, (long)d1, (long)d2); // R = P*Q+G;
    //     DP.write( DP.write(
    //         std::cout << "[(" , P ) << ") * (", Q) << ")]_" << d1 << '^' << d2 ;
    //     DP.write(std::cout << " = " , R) << std::endl;

    DP.axpy (T, P, Q, G);
    V=T;
    DP.truncin(V,(long)d1,(long)d2);

    S = G; U = G;
    DP.axpyin(S,P,Q,(long)d1,(long)d2);
    DP.axpyin(U,P,Q);
    W=U;
    DP.truncin(W,(long)d1,(long)d2);

    if( DP.areNEqual( V, R) || DP.areNEqual(W, S) || DP.areNEqual(V, W) ) {
        DP.write(std::cerr << "ERROR Axpy: ") << std::endl;
        std::cout << "  d1:= " << d1 << ';' << std::endl;
        std::cout << "  d2:= " << d2 << ';' << std::endl;
        DP.write(std::cout << "  P:= " , P) << ';' << std::endl;
        DP.write(std::cout << "  Q:= " , Q) << ';' << std::endl;
        DP.write(std::cout << "  G:= " , G) << ';' << std::endl;
        DP.write(std::cout << "  R:= " , R) << ';' << std::endl;
        DP.write(std::cout << "  T:= " , T) << ';' << std::endl;
        DP.write(std::cout << "  V:= " , V) << ';' << std::endl;
        DP.write(std::cout << "  S:= " , S) << ';' << std::endl;
        DP.write(std::cout << "  U:= " , U) << ';' << std::endl;
        DP.write(std::cout << "  W:= " , W) << ';' << std::endl;
        Degree dP; DP.degree(dP,P);
        Degree dQ; DP.degree(dQ,Q);
        Degree dG; DP.degree(dG,G);
        Degree vP; DP.val(vP,P);
        Degree vQ; DP.val(vQ,Q);
        Degree vG; DP.val(vG,G);
        std::cerr << "vP: " << vP << ", dP: " << dP << ", vQ: " << vQ << ", dQ: " << dQ << ", vG: " << vG << ", dG: " << dG << ", v: " << d1 << ", d: " << d2 << std::endl;
        return false;
    }
    return true;
}

bool TestAxmy(const TruncDom< GFqDom<int> >& DP, const TruncDom< GFqDom<int> >::Element& P, const TruncDom< GFqDom<int> >::Element& Q, const TruncDom< GFqDom<int> >::Element& G, size_t d1, size_t d2)
{
    ++TTcount;
    TruncDom< GFqDom<int> >::Element R, T, V, S, U, W;
    DP.axmy ( R, P, Q, G, (long)d1, (long)d2); // R = P*Q-G;
    //     DP.write( DP.write(
    //         std::cout << "[(" , P ) << ") + (", Q) << ")]_" << d1 << '^' << d2 ;
    //     DP.write(std::cout << " = " , R) << std::endl;

    DP.axmy (T, P, Q, G);
    V=T;
    DP.truncin(V,(long)d1,(long)d2);

    S = G; U = G;
    DP.axmyin(S,P,Q,(long)d1,(long)d2);
    DP.axmyin(U,P,Q);
    W=U;
    DP.truncin(W,(long)d1,(long)d2);

    if( DP.areNEqual( V, R) || DP.areNEqual(W, S) || DP.areNEqual(V, W) ) {
        std::cerr << "ERROR Axmy:" << std::endl;
        std::cout << "  d1:= " << d1 << ';' << std::endl;
        std::cout << "  d2:= " << d2 << ';' << std::endl;
        DP.write(std::cout << "  P:= " , P) << ';' << std::endl;
        DP.write(std::cout << "  Q:= " , Q) << ';' << std::endl;
        DP.write(std::cout << "  G:= " , G) << ';' << std::endl;
        DP.write(std::cout << "  R:= " , R) << ';' << std::endl;
        DP.write(std::cout << "  T:= " , T) << ';' << std::endl;
        DP.write(std::cout << "  V:= " , V) << ';' << std::endl;
        DP.write(std::cout << "  S:= " , S) << ';' << std::endl;
        DP.write(std::cout << "  U:= " , U) << ';' << std::endl;
        DP.write(std::cout << "  W:= " , W) << ';' << std::endl;
        Degree dP; DP.degree(dP,P);
        Degree dQ; DP.degree(dQ,Q);
        Degree vP; DP.val(vP,P);
        Degree vQ; DP.val(vQ,Q);
        std::cerr << "vR: " << vP << ", dR: " << dP << ", vP: " << vQ << ", dP: " << dQ << ", v: " << d1 << ", d: " << d2 << std::endl;
        return false;
    }
    return true;
}

bool TestMaxpy(const TruncDom< GFqDom<int> >& DP, const TruncDom< GFqDom<int> >::Element& P, const TruncDom< GFqDom<int> >::Element& Q, const TruncDom< GFqDom<int> >::Element& G, size_t d1, size_t d2)
{
    ++TTcount;
    TruncDom< GFqDom<int> >::Element R, T, V, S, U, W;
    DP.maxpy ( R, P, Q, G, (long)d1, (long)d2); // R = P*Q;
    //     DP.write( DP.write(
    //         std::cout << "[(" , P ) << ") + (", Q) << ")]_" << d1 << '^' << d2 ;
    //     DP.write(std::cout << " = " , R) << std::endl;

    DP.maxpy (T, P, Q, G);
    V=T;
    DP.truncin(V,(long)d1,(long)d2);

    S = G; U = G;
    DP.maxpyin(S,P,Q,(long)d1,(long)d2);
    DP.maxpyin(U,P,Q);
    W=U;
    DP.truncin(W,(long)d1,(long)d2);

    if( DP.areNEqual( V, R)  || DP.areNEqual(W, S) || DP.areNEqual(V, W) ) {
        std::cerr << "ERROR Maxpy:" << std::endl;
        DP.write(std::cout << "  R: " , R) << std::endl;
        DP.write(std::cout << "  T: " , T) << std::endl;
        DP.write(std::cout << "  V: " , V) << std::endl;
        Degree dP; DP.degree(dP,P);
        Degree dQ; DP.degree(dQ,Q);
        Degree vP; DP.val(vP,P);
        Degree vQ; DP.val(vQ,Q);
        std::cerr << "vR: " << vP << ", dR: " << dP << ", vP: " << vQ << ", dP: " << dQ << ", v: " << d1 << ", d: " << d2 << std::endl;
        return false;
    }
    return true;
}


int main(int argc, char ** argv) {

    long seed = (argc>1?atoi(argv[1]):BaseTimer::seed());
#ifdef __GIVARO_DEBUG
    std::cerr << "seed: " << seed << std::endl;
#endif

    GFqDom<int> Z101( 101, 1 );  // integers modulo 101

    // Polynomials over Z101, with X as indeterminate
    TruncDom< GFqDom<int> > DP101( Z101, Indeter("X") );
    TruncDom< GFqDom<int> >::Element P, Q, R, monomial;
    GFqDom<int>::Element tmp;

    DP101.assign( P, Z101.init(tmp,5) ); // P is degree 0 polynomial : 5 modulo 101
    DP101.init( monomial, Degree(1), 33U) ; // -33 X
    Degree deg,val;
    DP101.degree(deg,monomial);
    DP101.val(val,monomial);
    DP101.addin( P, monomial ); // P += monomial
    DP101.init( monomial, Degree(2), 12U) ; // 12 X^2
    DP101.addin( P, monomial ); // P is now 5-33*X+12*X^2

    //     // DP101.read( std::cin, P); // would read P as a succession of integers :
    //                                 // deg leadcoeff (lead-1)coeff ... unitcoeff
    Q = P;

    DP101.init( Q, Degree(0), 6U );
    DP101.init( monomial, Degree(4), 3U);
    DP101.addin( Q, monomial) ;
    DP101.init( monomial, Degree(1), 75U);
    DP101.addin( Q, monomial) ;
    DP101.init( monomial, Degree(3), 45U);
    DP101.subin( Q, monomial) ;
    // Q is now 3*X^4+75*X-45*X^3+6

    DP101.mulin( Q, Degree(15) ) ;
    DP101.mulin( monomial, Degree(32));
    DP101.addin( Q, monomial) ;
    DP101.subin( Q, monomial) ;

    DP101.divin( Q, Degree(5) ) ;
    DP101.divin( monomial, Degree(12));
    DP101.addin( Q, monomial) ;
    DP101.subin( Q, monomial) ;

    DP101.divin( monomial, Degree(10));
    DP101.addin( Q, monomial) ;
    DP101.subin( Q, monomial) ;

    DP101.divin( monomial, Degree(10));
    DP101.addin( Q, monomial) ;
    DP101.subin( Q, monomial) ;
    DP101.setval(Q);

    DP101.mulin( P, Degree(15) ) ;

    DP101.mul ( R, P, Q); // R = P*Q;
    DP101.mul ( R, P, Q, Degree(28), Degree(30)); // R = P*Q;
    DP101.mul ( R, P, Q, Degree(0), Degree(30)); // R = P*Q;
    DP101.mul ( R, P, Q, Degree(28), Degree(100)); // R = P*Q;
    DP101.mul ( R, P, Q, Degree(4), Degree(10)); // R = P*Q;
    DP101.mul ( R, P, Q, Degree(75), Degree(100)); // R = P*Q;
    DP101.setval(R);
    DP101.mulin( Q, Degree(3) ) ;

    DP101.add ( R, P, Q); // R = P*Q;




    bool success;

    success = TestAdd(DP101, P, Q, 28, 30);
    success &= TestAdd(DP101, P, Q, 0, 30);
    success &= TestAdd(DP101, P, Q, 28, 100);
    success &= TestAdd(DP101, P, Q, 4, 10);
    success &= TestAdd(DP101, P, Q, 75, 100);
    success &= TestAdd(DP101, P, Q, 12, 18);
    success &= TestAdd(DP101, P, Q, 13, 18);
    success &= TestAdd(DP101, P, Q, 14, 18);
    success &= TestAdd(DP101, P, Q, 15, 18);
    success &= TestAdd(DP101, P, Q, 16, 18);
    success &= TestAdd(DP101, P, Q, 17, 18);
    success &= TestAdd(DP101, P, Q, 18, 18);
    success &= TestAdd(DP101, P, Q, 12, 17);
    success &= TestAdd(DP101, P, Q, 13, 17);
    success &= TestAdd(DP101, P, Q, 14, 17);
    success &= TestAdd(DP101, P, Q, 15, 17);
    success &= TestAdd(DP101, P, Q, 16, 17);
    success &= TestAdd(DP101, P, Q, 17, 17);
    success &= TestAdd(DP101, P, Q, 12, 16);
    success &= TestAdd(DP101, P, Q, 13, 16);
    success &= TestAdd(DP101, P, Q, 14, 16);
    success &= TestAdd(DP101, P, Q, 15, 16);
    success &= TestAdd(DP101, P, Q, 16, 16);
    success &= TestAdd(DP101, P, Q, 12, 15);
    success &= TestAdd(DP101, P, Q, 13, 15);
    success &= TestAdd(DP101, P, Q, 14, 15);
    success &= TestAdd(DP101, P, Q, 15, 15);
    success &= TestAdd(DP101, P, Q, 12, 14);
    success &= TestAdd(DP101, P, Q, 13, 14);
    success &= TestAdd(DP101, P, Q, 14, 14);
    success &= TestAdd(DP101, P, Q, 12, 13);
    success &= TestAdd(DP101, P, Q, 13, 13);
    success &= TestAdd(DP101, P, Q, 12, 12);

    success &= TestAdd(DP101, Q, P, 12, 18);
    success &= TestAdd(DP101, Q, P, 13, 18);
    success &= TestAdd(DP101, Q, P, 14, 18);
    success &= TestAdd(DP101, Q, P, 15, 18);
    success &= TestAdd(DP101, Q, P, 16, 18);
    success &= TestAdd(DP101, Q, P, 17, 18);
    success &= TestAdd(DP101, Q, P, 18, 18);
    success &= TestAdd(DP101, Q, P, 12, 17);
    success &= TestAdd(DP101, Q, P, 13, 17);
    success &= TestAdd(DP101, Q, P, 14, 17);
    success &= TestAdd(DP101, Q, P, 15, 17);
    success &= TestAdd(DP101, Q, P, 16, 17);
    success &= TestAdd(DP101, Q, P, 17, 17);
    success &= TestAdd(DP101, Q, P, 12, 16);
    success &= TestAdd(DP101, Q, P, 13, 16);
    success &= TestAdd(DP101, Q, P, 14, 16);
    success &= TestAdd(DP101, Q, P, 15, 16);
    success &= TestAdd(DP101, Q, P, 16, 16);
    success &= TestAdd(DP101, Q, P, 12, 15);
    success &= TestAdd(DP101, Q, P, 13, 15);
    success &= TestAdd(DP101, Q, P, 14, 15);
    success &= TestAdd(DP101, Q, P, 15, 15);
    success &= TestAdd(DP101, Q, P, 12, 14);
    success &= TestAdd(DP101, Q, P, 13, 14);
    success &= TestAdd(DP101, Q, P, 14, 14);
    success &= TestAdd(DP101, Q, P, 12, 13);
    success &= TestAdd(DP101, Q, P, 13, 13);
    success &= TestAdd(DP101, Q, P, 12, 12);

    success &= TestAxpy(DP101, P,Q,R, 1, 2);
    success &= TestAxmy(DP101, P,Q,R, 1, 2);
    success &= TestMaxpy(DP101, P,Q,R, 1, 2);


    GivRandom generator((unsigned long)seed);

    for(size_t i=0; i<100; ++i) {
        long deg1 = (long)generator() % 75;
        long deg2 = (long)generator() % 85;
        long deg3 = (long)generator() % 155;
        long v1 = (long)generator() % 195;
        long v2 = v1 + (long)(generator() % 5);
        DP101.random(generator, P, Degree(deg1) );
        DP101.random(generator, Q, Degree(deg2) );
        DP101.random(generator, R, Degree(deg3) );
        Degree dP; DP101.degree(dP,P);
        Degree dQ; DP101.degree(dQ,Q);
        Degree dR; DP101.degree(dR,R);
        Degree vP; DP101.val(vP,P);
        Degree vQ; DP101.val(vQ,Q);
        Degree vR; DP101.val(vR,R);
        success &= TestAdd(DP101, P, Q, (size_t)v1, (size_t)v2 );
        success &= TestSub(DP101, P, Q, (size_t)v1, (size_t)v2 );
        success &= TestMul(DP101, P, Q, (size_t)v1, (size_t)v2 );
        success &= TestAxpy(DP101, P, Q, R, (size_t)v1, (size_t)v2 );
        success &= TestAxmy(DP101, P, Q, R, (size_t)v1, (size_t)v2 );
        success &= TestMaxpy(DP101, P, Q, R, (size_t)v1, (size_t)v2 );
    }

    success &= TestAxpy(DP101, P, Q, monomial, 11, 11);

#ifdef __GIVARO_DEBUG
    if (! success) {
        std::cerr << "Error: " << seed << std::endl;
    } else {
        std::cerr << "Success:" << TTcount << std::endl;
    }
#endif

    return (! success);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
