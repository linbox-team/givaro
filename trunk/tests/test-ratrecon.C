// Copyright(c)'1994-2009,2010 by The Givaro group
// This file is part of Givaro.
// BB: adapted,enhanced from examples/FiniteFields/ff-arith.C
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

#include <iostream>
#include <givaro/givpoly1.h>
#include <givaro/givfractiondomain.h>
#include <givaro/givzpz64std.h>

#ifdef GIVARO_DEBUG
long long TTcount = 0;
#endif

using namespace Givaro;

#define TEST_EQ( F, a, b ) \
if (!F.areEqual((a),(b))) {\
        F.write( F.write(std::cout,a) << "!=",b) << " failed (at line " <<  __LINE__ << ')' << std::endl; \
      return(-1); \
}

template<class PDomain>
int TestRR(PDomain& PolDom, GivRandom& generator, const Degree d) {

//     std::cout << "Start: " << d << std::endl;


    typename PDomain::Element P,Q,R,S,M,A,B, G;
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


    PolDom.invmod(R,Q,M);
// PolDom.write(std::cout << "R1:= ", R) << ';' << std::endl;
    PolDom.mulin(R, P);
// PolDom.write(std::cout << "R2:= ", R) << ';' << std::endl;
    PolDom.modin(R, M);
// PolDom.write(std::cout << "R3:= ", R) << ';' << std::endl;
    PolDom.ratrecon(A,B,R,M,d);

    typename PDomain::Type_t lc;
    PolDom.leadcoef(lc, B);
    PolDom.divin(A, lc);
    PolDom.divin(B, lc);

    PolDom.invmod(S,B,M);
// PolDom.write(std::cout << "S1:= ", S) << ';' << std::endl;
    PolDom.mulin(S, A);
// PolDom.write(std::cout << "S2:= ", S) << ';' << std::endl;
    PolDom.modin(S, M);
// PolDom.write(std::cout << "S3:= ", S) << ';' << std::endl;


//     PolDom.write(std::cout << "M:= ", M) << ';' << std::endl;
//     PolDom.write(std::cout << "P:= ", P) << ';' << std::endl;
//     PolDom.write(std::cout << "Q:= ", Q) << ';' << std::endl;
//     PolDom.write(std::cout << "R:= ", R) << ';' << std::endl;
//     PolDom.write(std::cout << "A:= ", A) << ';' << std::endl;
//     PolDom.write(std::cout << "B:= ", B) << ';' << std::endl;
//     PolDom.write( PolDom.write(std::cout << "P/Q: ", P) << "/", Q) << std::endl;
//     PolDom.write( PolDom.write(std::cout << "A/B: ", A) << "/", B) << std::endl;


//     TEST_EQ(PolDom, A, P);
//     TEST_EQ(PolDom, B, Q);
    TEST_EQ(PolDom, S, R);

#ifdef GIVARO_DEBUG
    ++TTcount;
#endif
    return 0;

}




int main(int argc, char ** argv)
{

    int seed = (argc>1?atoi(argv[1]):BaseTimer::seed());
#ifdef GIVARO_DEBUG
    std::cerr << "seed: " << seed << std::endl;
#endif
    Integer::seeding(seed);
    GivRandom generator(seed);

    typedef ZpzDom<Std64> Field;
    typedef Poly1Dom< Field, Dense > PolyZpz;
    typedef FracDom<PolyZpz> FracZpz;
    typedef Poly1Dom< FracZpz, Dense > PolyFracZpz;

    Field F101(101);
    Field F2(2);
    Field F65521(65521);
    bool success = true;

    {
        PolyZpz PZ(F101,"X");
        for(size_t j=0; j<4; ++j)
            for(Degree i=1; i<Degree(30); ++i)
                success &= (! TestRR(PZ, generator, i) );

        PolyFracZpz PFZ(PZ,"Y");
        for(size_t j=0; j<2; ++j)
            for(Degree i=2; i<Degree(5); ++i)
                success &= (! TestRR(PFZ, generator, i) );

    }

    {
        PolyZpz PZ(F2,"X");
        for(size_t j=0; j<4; ++j)
            for(Degree i=5; i<Degree(50); ++i)
                success &= (! TestRR(PZ, generator, i) );

        PolyFracZpz PFZ(PZ,"Y");
        for(size_t j=0; j<3; ++j)
            for(Degree i=7; i<Degree(15); ++i)
                success &= (! TestRR(PFZ, generator, i) );

    }

    {
        PolyZpz PZ(F65521,"X");
        for(size_t j=0; j<5; ++j)
            for(Degree i=1; i<Degree(30); ++i)
                success &= (! TestRR(PZ, generator, i) );
        PolyFracZpz PFZ(PZ,"Y");
        for(size_t j=0; j<2; ++j)
            for(Degree i=2; i<Degree(5); ++i)
                success &= (! TestRR(PFZ, generator, i) );

    }





#ifdef GIVARO_DEBUG
    if (! success) {
        std::cerr << "Error: " << seed << std::endl;
    } else {
        std::cerr << "Success:" << TTcount << std::endl;
    }
#endif

    return (! success);
}
