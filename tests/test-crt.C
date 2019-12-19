// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

/*! @file tests/test-crt.C
 * @ingroup tests
 * @brief NO DOC
 */

#include <iostream>
#include <givaro/givintprime.h>
#include <givaro/montgomery.h>
#include <givaro/modular.h>
#include <givaro/gfq.h>
#include <givaro/chineseremainder.h>    // Chinese Remainder of two elements
#include <givaro/givrns.h>    // Chinese Remainder of an array of elements
#include <givaro/givrnsfixed.h>    // Chinese Remainder with fixed primes
#include <givaro/givtimer.h>
#include <givaro/givrandom.h>

using namespace Givaro;

typedef GFqDom<int32_t>            Field1;

typedef Modular<int16_t>        Field2;
typedef Modular<int32_t>        Field3;
typedef Modular<int64_t>        Field4;
typedef Modular<uint16_t>       Field5;
typedef Modular<uint32_t>       Field6;
typedef Modular<uint64_t>       Field7;
typedef Modular<Log16>          Field8;
typedef Modular<Integer>        Field9;

typedef Montgomery<int32_t>     Field10;

//typedef Modular<int8_t>         Field11; Too small for the test
//typedef Modular<uint8_t>        Field12; Too small for the test

template <typename Field>
Integer tmain(int argc, char ** argv, const GivRandom& generator, const Integer UpperBoundChard, bool isFieldModular=true) {


    typedef RNSsystem<Integer, Field>      CRTSystem;
    typedef typename CRTSystem::domains      Domains;
    typedef typename CRTSystem::array     Elements;

    typedef RNSsystemFixed<Integer>    CRTSystemFixed;
    typedef CRTSystemFixed::array             Prime_t;

    IntPrimeDom ID;
    Prime_t  Primes( argc>1 ? (size_t)atoi(argv[1]):15);
    Elements Moduli( Primes.size() );
    Prime_t  ModuliInts( Primes.size() );
    Domains PrimeDoms( Primes.size() );
    auto p = Primes.begin();
    auto i = PrimeDoms.begin();
    auto e = Moduli.begin();
    auto m = ModuliInts.begin();
    Integer b, M(1), a(UpperBoundChard);

    {

        for(; i != PrimeDoms.end(); ++i, ++e, ++p, ++m) {
            *i = Field( ID.prevprimein( a ) );
            *p = a;

            assert(i->characteristic() == *p);

            uint64_t tmp = generator();
            i->init(*e,  tmp );
            i->convert(*m,*e);
            M *= a;
        }
    }

    CRTSystem CRT( PrimeDoms );
    assert(CRT.size() == Primes.size());

    Timer tim; tim.clear(); tim.start();
    Integer TestRes, tmp; CRT.RnsToRing( TestRes, Moduli );
    tim.stop();


#ifdef __GIVARO_DEBUG
    std::cerr << "MAXC: " << Field::maxCardinality() << std::endl;
    Field(PrimeDoms.front()).write( std::cerr << tim << " using ") << std::endl;
#endif

#ifdef __GIVARO_DEBUG
    if (PrimeDoms.size() < 50) {
        i = PrimeDoms.begin();
        e = Moduli.begin();
        for( ; i != PrimeDoms.end(); ++i, ++e)
            i->write(std::cout << TestRes << " mod " << i->characteristic() << " = ", *e) << ";" << std::endl;
    }
#endif

    if (isFieldModular) {
        Timer timf;
        timf.clear(); timf.start();
        CRTSystemFixed CRTFixed( Primes );
        timf.stop();
#ifdef __GIVARO_DEBUG
        std::cerr << "CRTFixed init : " << timf << std::endl;
#endif

        timf.clear(); timf.start();
        CRTFixed.RnsToRing( b, ModuliInts );
        timf.stop();
#ifdef __GIVARO_DEBUG
        std::cerr << "CRTFixed : " << timf << std::endl;
#endif

        if (TestRes != b) {
            std::cerr << "Error Field: ";
            PrimeDoms[0].write(std::cerr) << std::endl;
            std::cerr << "incoherency between normal : " << TestRes
            << " and fixed : " << b << std::endl;
            exit(1); // très peu probable que res=0 pour de vrai.
        }

    }


    Elements Verifs( Primes.size() );
    CRT.RingToRns( Verifs, TestRes );

    // #ifdef __GIVARO_DEBUG
    //     for (const auto v : Verifs)
    // 	std::cerr << v << std::endl;
    // #endif

    typename Elements::const_iterator v = Verifs.begin();
    i = PrimeDoms.begin();
    e = Moduli.begin();
    for( ; i != PrimeDoms.end(); ++i, ++e, ++v) {
        if (! i->areEqual(*e, *v) ) {
            i->write( i->write( std::cerr << "Error e: ", *e) << " != ", *v) << std::endl;
            i->write( std::cerr << "incoherency within ") << std::endl;
            exit(2);
        }
    }

    Integer pr( generator() >>(argc>2?atoi(argv[2]):17) ), res;
    if ( (Field::maxCardinality() > 0) && (pr > Field::maxCardinality() ) ) {
        pr = Field::maxCardinality()/2;
    }
    Field F( ID.prevprimein(pr) );
    typename Field::Element el;
    F.init(el, generator() );

    ChineseRemainder<IntPrimeDom, Field> CRA(ID, M, F);
    CRA( res, TestRes, el);

    ID.mod(tmp,res,M);
    if (! ID.areEqual(tmp,TestRes)) {
        std::cerr << "Error CRA: " << res << " mod " << M << " != " << TestRes << ";"  << std::endl;
    }


#ifdef __GIVARO_DEBUG
    std::cout << res << " mod " << M << " = " << TestRes << ";"  << std::endl;
    std::cout << res << " mod " << F.characteristic() << " = " << F.convert(tmp, el) << ";"  << std::endl;
#endif

    return TestRes;
}



int main(int argc, char ** argv)
{
    ::Givaro::GivaroMain::Init();
#ifdef __GIVARO_DEBUG
    Givaro::GivMMInfo MemoryInfo;
#endif

    // argv[1] : number of primes
    // argv[2] : 2^{32-j} is size of primes
    // argv[3] : seed for generator

    GivRandom seedor( argc>3 ? (unsigned)atoi(argv[3]): (unsigned)BaseTimer::seed() );
    uint64_t seed = seedor.seed();
    GivRandom generator(seed);

    Integer ubc( generator() >>(argc>2?atoi(argv[2]):17) );
    if ( (Field4::maxCardinality() > 0)
         && ( ubc > Field4::maxCardinality() ) ) {
        ubc = Field4::maxCardinality();
    }
    if (ubc<59) ubc = 59; // at least 15 primes
    Integer a4 = tmain<Field4>(argc, argv, GivRandom(seed), ubc);
    Integer a7 = tmain<Field7>(argc, argv, GivRandom(seed), ubc);
    Integer a9 = tmain<Field9>(argc, argv, GivRandom(seed), ubc);

    if ( ubc > Field3::maxCardinality() )  {
        ubc = Field3::maxCardinality();
    }
    if (ubc<59) ubc = 59; // at least 15 primes

    Integer a1 = tmain<Field1>(argc, argv, GivRandom(seed), ubc, false);
    Integer a3 = tmain<Field3>(argc, argv, GivRandom(seed), ubc);
    Integer a6 = tmain<Field6>(argc, argv, GivRandom(seed), ubc);

    if ( ubc > Field8::maxCardinality() )  {
        ubc = Field8::maxCardinality();
    }
    if (ubc<59) ubc = 59; // at least 15 primes

    Integer a8 = tmain<Field8>(argc, argv, GivRandom(seed), ubc);
    Integer a10 = tmain<Field10>(argc, argv, GivRandom(seed), ubc, false);


    if ( ubc > Field2::maxCardinality() )  {
        ubc = Field2::maxCardinality();
    }
    if (ubc<59) ubc = 59; // at least 15 primes

    Integer a2 = tmain<Field2>(argc, argv, GivRandom(seed), ubc);
    Integer a5 = tmain<Field5>(argc, argv, GivRandom(seed), ubc);

#ifdef __GIVARO_DEBUG
    std::cerr << "seed: " << seed << std::endl;
    std::cerr << "a1: " << a1 << std::endl;
    std::cerr << "a3: " << a3 << std::endl;
    std::cerr << "a4: " << a4 << std::endl;
    std::cerr << "a6: " << a6 << std::endl;
    std::cerr << "a7: " << a7 << std::endl;
    std::cerr << "a8: " << a8 << std::endl;
    std::cerr << "a9: " << a9 << std::endl;
    std::cerr << "a10: " << a10 << std::endl;
    std::cerr << "a2: " << a2 << std::endl;
    std::cerr << "a5: " << a5 << std::endl;
#endif

    bool success = true;
    bool suctest = true;
    success = (a4 == a7); suctest &= success;
    if (! success) std::cerr << "ERROR a4 != a7" << std::endl;
    success = (a4 == a9); suctest &= success;
    if (! success) std::cerr << "ERROR a4 != a9" << std::endl;

    success = (a1 == a3); suctest &= success;
    if (! success) std::cerr << "ERROR a1 != a3" << std::endl;
    success = (a1 == a6); suctest &= success;
    if (! success) std::cerr << "ERROR a1 != a6" << std::endl;

    success = (a8 == a10); suctest &= success;
    if (! success) std::cerr << "ERROR a8 != a10" << std::endl;

    success = (a2 == a5); suctest &= success;
    if (! success) std::cerr << "ERROR a2 != a5" << std::endl;

#ifdef __GIVARO_DEBUG
    if (! suctest)
        std::cerr << "Error: " << seed << std::endl;
    MemoryInfo.print(std::cerr) << std::endl;
#endif

    ::Givaro::GivaroMain::End();
    return (! suctest);
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
