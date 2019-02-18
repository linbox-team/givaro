// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

/*! @file examples/Polynomial/trunc_arith.C
 * @ingroup examples
 * @ingroup polynomials
 * @example examples/Polynomial/trunc_arith.C
 * @brief NO DOC
 */
#include <iostream>
#include <givaro/givrandom.h>
#include <givaro/givtimer.h>
#include <givaro/gfq.h>
#include <givaro/givpoly1.h>
#include <givaro/givtruncdomain.h>

using namespace Givaro;



long long TTcount = 0;

bool TestAdd(const TruncDom< GFqDom<int> >& DP, const TruncDom< GFqDom<int> >::Element& P, const TruncDom< GFqDom<int> >::Element& Q, size_t d1, size_t d2)
{
    ++TTcount;
    TruncDom< GFqDom<int> >::Element R, T, V;
    DP.add ( R, P, Q, (int64_t)d1, (int64_t)d2); // R = P*Q;
    //     DP.write( DP.write(
    //         std::cout << "[(" , P ) << ") + (", Q) << ")]_" << d1 << '^' << d2 ;
    //     DP.write(std::cout << " = " , R) << std::endl;

    DP.add (T, P, Q);
    V=T;
    DP.truncin(V,(int64_t)d1,(int64_t)d2);

    if( DP.areNEqual( V, R) ) {
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
    TruncDom< GFqDom<int> >::Element R, T, V;
    DP.sub ( R, P, Q, (int64_t)d1, (int64_t)d2); // R = P*Q;
    //     DP.write( DP.write(
    //         std::cout << "[(" , P ) << ") + (", Q) << ")]_" << d1 << '^' << d2 ;
    //     DP.write(std::cout << " = " , R) << std::endl;

    DP.sub (T, P, Q);
    V=T;
    DP.truncin(V,(int64_t)d1,(int64_t)d2);

    if( DP.areNEqual( V, R) ) {
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
    TruncDom< GFqDom<int> >::Element R, T, V;
    DP.mul ( R, P, Q, (int64_t)d1, (int64_t)d2); // R = P*Q;
    //     DP.write( DP.write(
    //         std::cout << "[(" , P ) << ") + (", Q) << ")]_" << d1 << '^' << d2 ;
    //     DP.write(std::cout << " = " , R) << std::endl;

    DP.mul (T, P, Q);
    V=T;
    DP.truncin(V,(int64_t)d1,(int64_t)d2);

    if( DP.areNEqual( V, R) ) {
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
    TruncDom< GFqDom<int> >::Element R, T, V;
    DP.axpy ( R, P, Q, G, (int64_t)d1, (int64_t)d2); // R = P*Q;
    //     DP.write( DP.write(
    //         std::cout << "[(" , P ) << ") + (", Q) << ")]_" << d1 << '^' << d2 ;
    //     DP.write(std::cout << " = " , R) << std::endl;

    DP.axpy (T, P, Q, G);
    V=T;
    DP.truncin(V,(int64_t)d1,(int64_t)d2);

    if( DP.areNEqual( V, R) ) {
        std::cerr << "ERROR Axpy:" << std::endl;
        DP.write(std::cout << "  P: " , P) << std::endl;
        DP.write(std::cout << "  Q: " , Q) << std::endl;
        DP.write(std::cout << "  G: " , G) << std::endl;
        DP.write(std::cout << "  R: " , R) << std::endl;
        DP.write(std::cout << "  T: " , T) << std::endl;
        DP.write(std::cout << "  V: " , V) << std::endl;
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
    TruncDom< GFqDom<int> >::Element R, T, V;
    DP.axmy ( R, P, Q, G, (int64_t)d1, (int64_t)d2); // R = P*Q;
    //     DP.write( DP.write(
    //         std::cout << "[(" , P ) << ") + (", Q) << ")]_" << d1 << '^' << d2 ;
    //     DP.write(std::cout << " = " , R) << std::endl;

    DP.axmy (T, P, Q, G);
    V=T;
    DP.truncin(V,(int64_t)d1,(int64_t)d2);

    if( DP.areNEqual( V, R) ) {
        std::cerr << "ERROR Axmy:" << std::endl;
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

bool TestMaxpy(const TruncDom< GFqDom<int> >& DP, const TruncDom< GFqDom<int> >::Element& P, const TruncDom< GFqDom<int> >::Element& Q, const TruncDom< GFqDom<int> >::Element& G, size_t d1, size_t d2)
{
    ++TTcount;
    TruncDom< GFqDom<int> >::Element R, T, V;
    DP.maxpy ( R, P, Q, G, (int64_t)d1, (int64_t)d2); // R = P*Q;
    //     DP.write( DP.write(
    //         std::cout << "[(" , P ) << ") + (", Q) << ")]_" << d1 << '^' << d2 ;
    //     DP.write(std::cout << " = " , R) << std::endl;

    DP.maxpy (T, P, Q, G);
    V=T;
    DP.truncin(V,(int64_t)d1,(int64_t)d2);

    if( DP.areNEqual( V, R) ) {
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

    {
        int64_t seed = (argc>1?atoi(argv[1]):BaseTimer::seed());
        std::cerr << "seed: " << seed << std::endl;

        GFqDom<int> Z101( 101, 1 );  // integers modulo 101

        // Polynomials over Z101, with X as indeterminate
        TruncDom< GFqDom<int> > DP101( Z101, Indeter("X") );
        TruncDom< GFqDom<int> >::Element P, Q, R, monomial;
        GFqDom<int>::Element tmp;

        DP101.assign( P, Z101.init(tmp,5) ); // P is degree 0 polynomial : 5 modulo 101
        DP101.write( std::cout << "P: " , P )<< std::endl;
        DP101.init( monomial, Degree(1), 33U) ; // -33 X
        DP101.write( std::cout << "m: " , monomial )<< std::endl;
        Degree deg,val;
        DP101.degree(deg,monomial);
        DP101.val(val,monomial);
        DP101.write( std::cout << "[m]_" << val << '^' << deg << ": " , monomial )<< std::endl;
        DP101.addin( P, monomial ); // P += monomial
        DP101.write( std::cout << "P: " , P )<< std::endl;
        DP101.init( monomial, Degree(2), 12U) ; // 12 X^2
        DP101.write( std::cout << "m: " , monomial )<< std::endl;
        DP101.addin( P, monomial ); // P is now 5-33*X+12*X^2
        DP101.write( std::cout << "P: " , P )<< std::endl;

        //     // DP101.read( std::cin, P); // would read P as a succession of integers :
        //                                 // deg leadcoeff (lead-1)coeff ... unitcoeff
        Q = P;
        DP101.write( std::cout << "P: " , P )<< std::endl;
        DP101.write( std::cout << "Q: " , Q )<< std::endl;


        DP101.init( Q, Degree(0), 6U );
        DP101.write( std::cout << "Q: " , Q )<< std::endl;
        DP101.init( monomial, Degree(4), 3U);
        DP101.write( std::cout << "m: " , monomial )<< std::endl;
        DP101.addin( Q, monomial) ;
        DP101.write( std::cout << "Q: " , Q )<< std::endl;
        DP101.init( monomial, Degree(1), 75U);
        DP101.write( std::cout << "m: " , monomial )<< std::endl;
        DP101.addin( Q, monomial) ;
        DP101.write( std::cout << "Q: " , Q )<< std::endl;
        DP101.init( monomial, Degree(3), 45U);
        DP101.write( std::cout << "m: " , monomial )<< std::endl;
        DP101.subin( Q, monomial) ;
        DP101.write( std::cout << "Q: " , Q )<< std::endl;
        // Q is now 3*X^4+75*X-45*X^3+6

        DP101.mulin( Q, Degree(15) ) ;
        DP101.write( std::cout << "Q: " , Q )<< std::endl;
        DP101.mulin( monomial, Degree(32));
        DP101.write( std::cout << "m: " , monomial )<< std::endl;
        DP101.addin( Q, monomial) ;
        DP101.write( std::cout << "Q: " , Q )<< std::endl;
        DP101.subin( Q, monomial) ;
        DP101.write( std::cout << "Q: " , Q )<< std::endl;

        DP101.divin( Q, Degree(5) ) ;
        DP101.write( std::cout << "Q: " , Q )<< std::endl;
        DP101.divin( monomial, Degree(12));
        DP101.write( std::cout << "m: " , monomial )<< std::endl;
        DP101.addin( Q, monomial) ;
        DP101.write( std::cout << "Q: " , Q )<< std::endl;
        DP101.subin( Q, monomial) ;
        DP101.write( std::cout << "Q: " , Q )<< std::endl;

        DP101.divin( monomial, Degree(10));
        DP101.write( std::cout << "m: " , monomial )<< std::endl;
        DP101.addin( Q, monomial) ;
        DP101.write( std::cout << "Q: " , Q )<< std::endl;
        DP101.subin( Q, monomial) ;
        DP101.write( std::cout << "Q: " , Q )<< std::endl;

        DP101.divin( monomial, Degree(10));
        DP101.write( std::cout << "m: " , monomial )<< std::endl;
        DP101.addin( Q, monomial) ;
        DP101.write( std::cout << "Q: " , Q )<< std::endl;
        DP101.subin( Q, monomial) ;
        DP101.write( std::cout << "Q: " , Q )<< std::endl;
        DP101.setval(Q);
        DP101.write( std::cout << "Q: " , Q )<< std::endl;

        DP101.mulin( P, Degree(15) ) ;

        DP101.mul ( R, P, Q); // R = P*Q;
        DP101.write( DP101.write(
                                 std::cout << "(" , P ) << ") * (", Q) << ")";
        DP101.write(std::cout << " = " , R) << std::endl;

        DP101.mul ( R, P, Q, Degree(28), Degree(30)); // R = P*Q;
        DP101.write( DP101.write(
                                 std::cout << "[(" , P ) << ") * (", Q) << ")]_28^30";
        DP101.write(std::cout << " = " , R) << std::endl;
        DP101.mul ( R, P, Q, Degree(0), Degree(30)); // R = P*Q;
        DP101.write( DP101.write(
                                 std::cout << "[(" , P ) << ") * (", Q) << ")]_0^30";
        DP101.write(std::cout << " = " , R) << std::endl;
        DP101.mul ( R, P, Q, Degree(28), Degree(100)); // R = P*Q;
        DP101.write( DP101.write(
                                 std::cout << "[(" , P ) << ") * (", Q) << ")]_28^100";
        DP101.write(std::cout << " = " , R) << std::endl;
        DP101.mul ( R, P, Q, Degree(4), Degree(10)); // R = P*Q;
        DP101.write( DP101.write(
                                 std::cout << "[(" , P ) << ") * (", Q) << ")]_4^10";
        DP101.write(std::cout << " = " , R) << std::endl;
        DP101.mul ( R, P, Q, Degree(75), Degree(100)); // R = P*Q;
        DP101.write( DP101.write(
                                 std::cout << "[(" , P ) << ") * (", Q) << ")]_75^100";
        DP101.write(std::cout << " = " , R) << std::endl;
        DP101.setval(R);
        DP101.write( std::cout << "R: " , R )<< std::endl;

        DP101.mulin( Q, Degree(3) ) ;

        DP101.add ( R, P, Q); // R = P*Q;
        DP101.write( DP101.write(
                                 std::cout << "(" , P ) << ") + (", Q) << ")";
        DP101.write(std::cout << " = " , R) << std::endl;

        /*
           DP101.gcd ( R, P, Q); //

           DP101.write( DP101.write( DP101.write(
           std::cout << "gcd(", P ) << ",", Q) << ") = ", R) << std::endl;

           DP101.lcm ( R, P, Q); //
           DP101.write( DP101.write( DP101.write(
           std::cout << "lcm(", P ) << ",", Q) << ") = ", R) << std::endl;
           DP101.lcm ( R, Q, P); //
           DP101.write( DP101.write( DP101.write(
           std::cout << "lcm(", Q ) << ",", P) << ") = ", R) << std::endl;

*/
    }
    return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
