// Copyright(c)'1994-2009,2010 by The Givaro group
// This file is part of Givaro.
// BB: adapted,enhanced from examples/FiniteFields/ff-arith.C
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

#include <iostream>
#include <givaro/givinteger.h>
#include <givaro/givpoly1.h>
#include <givaro/givfractiondomain.h>
#include <givaro/modular-integral.h>

#ifdef __GIVARO_DEBUG
long long TTcount = 0;
#endif

using namespace Givaro;

#define TEST_EQ( F, a, b ) \
if (!F.areEqual((a),(b))) {\
    F.write( F.write(std::cout,a) << "!=",b) << " failed (at line " <<  __LINE__ << ')' << std::endl; \
    return(-1); \
}




template<class RingDomain, class BoundingStathme>
int TestRR(RingDomain& RDom, GivRandom& generator,
           const typename RingDomain::Element& M,
           const typename RingDomain::Element& P,
           const typename RingDomain::Element& Q,
           const BoundingStathme& d) {
//     RDom.write(std::clog << "M:= ", M) << ';' << std::endl;
//     RDom.write(std::clog << "P:= ", P) << ';' << std::endl;
//     RDom.write(std::clog << "Q:= ", Q) << ';' << std::endl;
//     std::clog << "Bound:=" << d  << ';' << std::endl;

    typename RingDomain::Element R,S,A,B;

    RDom.invmod(R,Q,M);
//     RDom.write(std::clog << "R1:= ", R) << ';' << std::endl;
    RDom.mulin(R, P);
//     RDom.write(std::clog << "R2:= ", R) << ';' << std::endl;
    RDom.modin(R, M);
//     RDom.write(std::clog << "R3:= ", R) << ';' << std::endl;
    RDom.ratrecon(A,B,R,M,d,true);

//     RDom.write(std::clog << "A:= ", A) << ';' << std::endl;
//     RDom.write(std::clog << "B:= ", B) << ';' << std::endl;

    RDom.mul(S, R, B);
//     RDom.write(std::clog << "S1:= ", S) << ';' << std::endl;
    RDom.modin(S, M);
//     RDom.write(std::clog << "S2:= ", S) << ';' << std::endl;
    RDom.subin(S, A);
//     RDom.write(std::clog << "S3:= ", S) << ';' << std::endl;
    RDom.modin(S, M);
//     RDom.write(std::clog << "S4:= ", S) << ';' << std::endl;

//     RDom.write(std::clog << "M:= ", M) << ';' << std::endl;
//     RDom.write(std::clog << "P:= ", P) << ';' << std::endl;
//     RDom.write(std::clog << "Q:= ", Q) << ';' << std::endl;
//     RDom.write(std::clog << "R:= ", R) << ';' << std::endl;
//     RDom.write(std::clog << "A:= ", A) << ';' << std::endl;
//     RDom.write(std::clog << "B:= ", B) << ';' << std::endl;
//     RDom.write(std::clog << "S:= ", S) << ';' << std::endl;
//     RDom.write( RDom.write(std::clog << "P/Q: ", P) << "/", Q) << std::endl;
//     RDom.write( RDom.write(std::clog << "A/B: ", A) << "/", B) << std::endl;

    //     TEST_EQ(RDom, A, P);
    //     TEST_EQ(RDom, B, Q);
    TEST_EQ(RDom, S, RDom.zero);

#ifdef __GIVARO_DEBUG
    ++TTcount;
#endif
    return 0;

}





int TestRatRR(GivRandom& generator, const size_t b) {

//     std::clog << "Start: " << b << std::endl;


    ZRing<Integer> IntDom;

    Integer M, P, Q, G;
    IntDom.random(generator, M, (b<<1) + 2 );

    do {
        IntDom.random(generator, Q, b);
        IntDom.gcd(G, M, Q);
    } while( G > 1);
    do {
        IntDom.random(generator, P, b-1);
        IntDom.gcd(G, P, Q);
    } while( G > 1);

    return TestRR(IntDom, generator, M, P, Q, P+1);
}



template<class PDomain>
int TestPolRR(PDomain& PolDom, GivRandom& generator, const Degree d) {

//     std::clog << "Start: " << d << std::endl;


    typename PDomain::Element M, P, Q, G;
    PolDom.random(generator, M, (d << 1) + 2);

    Degree dG;
    do {
        PolDom.random(generator, Q, d);
        PolDom.getdomain().assign(Q[Q.size()-1], PolDom.getdomain().one);
        PolDom.gcd(G, M, Q);
        PolDom.degree(dG, G);
    } while( dG > 0);
    do {
        PolDom.random(generator, P, d-1);
        PolDom.gcd(G, P, Q);
        PolDom.degree(dG, G);
    } while( dG > 0);


    return TestRR(PolDom, generator, M, P, Q, d);
}


int main(int argc, char ** argv)
{

    int seed = int(argc>1?atoi(argv[1]):BaseTimer::seed());
#ifdef __GIVARO_DEBUG
    std::cerr << "seed: " << seed << std::endl;
#endif
    Integer::seeding((unsigned long)seed);
    GivRandom generator((unsigned long)seed);


//     RDom.write(std::clog << "R3:= ", R) << ';' << std::endl;
    ZRing<Integer> IntDom;
    Integer R(75), M(250), A, B, d1(17), d2(26);

    bool success = true;

        // Should be impossible;
    bool reconstructed = IntDom.ratrecon(A,B,R,M,d1);
    success &= (!reconstructed);

        // Now it should be ok
    reconstructed = IntDom.ratrecon(A,B,R,M,d2);
    success &= (reconstructed && (IntDom.isZero((R*B-A)%M)));

    reconstructed = IntDom.ratrecon(A,B,R,M,d1,false);
    success &= (reconstructed && (IntDom.isZero((R*B-A)%M)));


    typedef Modular<int64_t> Field;
    typedef Poly1Dom< Field, Dense > PolyZpz;
    typedef FracDom<PolyZpz> FracZpz;
    typedef Poly1Dom< FracZpz, Dense > PolyFracZpz;

    Field F101(101);
    Field F2(2);
    Field F65521(65521);

    {
        for(size_t loop=0; loop<100; ++loop)
            for (size_t b=10; b<10000; b <<=1)
                success &= (! TestRatRR(generator, b) );
    }



    {
        PolyZpz PZ(F101,"X");
        for(size_t j=0; j<4; ++j)
            for(Degree i=1; i<Degree(30); ++i)
                success &= (! TestPolRR(PZ, generator, i) );

        PolyFracZpz PFZ(PZ,"Y");
        for(size_t j=0; j<2; ++j)
            for(Degree i=2; i<Degree(5); ++i)
                success &= (! TestPolRR(PFZ, generator, i) );

    }

    {
        PolyZpz PZ(F2,"X");
        for(size_t j=0; j<4; ++j)
            for(Degree i=5; i<Degree(50); ++i)
                success &= (! TestPolRR(PZ, generator, i) );

        PolyFracZpz PFZ(PZ,"Y");
        for(size_t j=0; j<3; ++j)
            for(Degree i=7; i<Degree(15); ++i)
                success &= (! TestPolRR(PFZ, generator, i) );

    }

    {
        PolyZpz PZ(F65521,"X");
        for(size_t j=0; j<5; ++j)
            for(Degree i=1; i<Degree(30); ++i)
                success &= (! TestPolRR(PZ, generator, i) );
        PolyFracZpz PFZ(PZ,"Y");
        for(size_t j=0; j<2; ++j)
            for(Degree i=2; i<Degree(5); ++i)
                success &= (! TestPolRR(PFZ, generator, i) );

    }



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
