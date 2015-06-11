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

typedef GFqDom<long>            Field1;

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
Integer tmain(int argc, char ** argv, const GivRandom& generator)
{
    typedef RNSsystem<Integer, Field>      CRTSystem;
    typedef typename CRTSystem::domains      Domains;
    typedef typename CRTSystem::array     Elements;

    typedef RNSsystemFixed<Integer>    CRTSystemFixed;
    typedef CRTSystemFixed::array             Prime_t;

    IntPrimeDom ID;
    Integer a( generator() >>(argc>2?atoi(argv[2]):17) ), M(1), b;

    Prime_t  Primes( argc>1 ? (size_t)atoi(argv[1]):15);
    Elements Moduli( Primes.size() );
    Prime_t  ModuliInts( Primes.size() );
    Domains PrimeDoms( Primes.size() );

    auto p = Primes.begin();
    auto i = PrimeDoms.begin();
    auto e = Moduli.begin();
    auto m = ModuliInts.begin();
    
    for(; i != PrimeDoms.end(); ++i, ++e, ++p, ++m) {
        *i = Field( ID.nextprimein( a ) );
        *p = a;
	assert(i->characteristic() == *p);
        i->init(*e,  generator() );
        i->convert(*m,*e);
        M *= a;
    }

    CRTSystem CRT( PrimeDoms );
    assert(CRT.size() == Primes.size());

    Timer tim; tim.clear(); tim.start();
    CRT.RnsToRing( a, Moduli );
    tim.stop();
#ifdef GIVARO_DEBUG
    Field(PrimeDoms.front()).write( std::cerr << tim << " using ") << std::endl;
#endif

#ifdef GIVARO_DEBUG
    if (PrimeDoms.size() < 50) {
        i = PrimeDoms.begin();
        e = Moduli.begin();
        for( ; i != PrimeDoms.end(); ++i, ++e)
            i->write(std::cout << a << " mod " << i->characteristic() << " = ", *e) << ";" << std::endl;
    }
#endif

    Timer timf;
    timf.clear(); timf.start();
    CRTSystemFixed CRTFixed( Primes );
    timf.stop();
#ifdef GIVARO_DEBUG
    std::cerr << "CRTFixed init : " << timf << std::endl;
#endif

    timf.clear(); timf.start();
    CRTFixed.RnsToRing( b, ModuliInts );
    timf.stop();
#ifdef GIVARO_DEBUG
    std::cerr << "CRTFixed : " << timf << std::endl;
#endif

    if (a != b) {
        std::cerr << "Field: ";
        PrimeDoms[0].write(std::cerr) << std::endl;
        std::cerr << "incoherency between normal : " << a
        << " and fixed : " << b << std::endl;
        return 0 ; // très peu probable que res=0 pour de vrai.
    }

    Elements Verifs( Primes.size() );
    CRT.RingToRns( Verifs, a );

#ifdef GIVARO_DEBUG
    for (const auto v : Verifs)
	std::cerr << v << std::endl;
#endif

    typename Elements::const_iterator v = Verifs.begin();
    i = PrimeDoms.begin();
    e = Moduli.begin();
    for( ; i != PrimeDoms.end(); ++i, ++e, ++v) {
        if (! i->areEqual(*e, *v) ) {
            i->write( std::cerr << "incoherency within ") << std::endl;
            exit(1);
        }
    }

    Integer pr( generator() >>(argc>2?atoi(argv[2]):17) ), res;
    Field F( ID.nextprimein(pr) );
    typename Field::Element el;
    F.init(el, generator() );

    ChineseRemainder<IntPrimeDom, Field> CRA(ID, M, F);
    CRA( res, a, el);

#ifdef GIVARO_DEBUG
    std::cout << res << " mod " << M << " = " << a << ";"  << std::endl;
    std::cout << res << " mod " << F.characteristic() << " = " << F.convert(a, el) << ";"  << std::endl;
#endif

    return  res;
}



int main(int argc, char ** argv)
{
    ::Givaro::GivaroMain::Init();
#ifdef GIVARO_DEBUG
    Givaro::GivMMInfo MemoryInfo;
#endif

    // argv[1] : number of primes
    // argv[2] : 2^{32-j} is size of primes
    // argv[3] : seed for generator

    GivRandom seedor( argc>3 ? (unsigned)atoi(argv[3]): (unsigned)BaseTimer::seed() );
    unsigned long seed = seedor.seed();

    Integer a1 = tmain<Field1>(argc, argv, GivRandom(seed));
    Integer a2 = tmain<Field2>(argc, argv, GivRandom(seed));
    Integer a3 = tmain<Field3>(argc, argv, GivRandom(seed));
    Integer a4 = tmain<Field4>(argc, argv, GivRandom(seed));
    
    Integer a5 = tmain<Field5>(argc, argv, GivRandom(seed));
    Integer a6 = tmain<Field6>(argc, argv, GivRandom(seed));
    Integer a7 = tmain<Field7>(argc, argv, GivRandom(seed));
    Integer a8 = tmain<Field8>(argc, argv, GivRandom(seed));
    
    Integer a9 = tmain<Field9>(argc, argv, GivRandom(seed));
    Integer a10 = tmain<Field10>(argc, argv, GivRandom(seed));

    if (!(a1 & a2 & a3 & a4 & a5 & a6 & a7 & a8 & a9 & a10)) {
#ifdef GIVARO_DEBUG
        std::cerr << "one test failed" << std::endl;
#endif
        return false ;
    }

    bool success = true;
    success &= (a1 == a2);
    if (! success) std::cerr << "ERROR a1 != a2" << std::endl;
    success &= (a3 == a4);
    if (! success) std::cerr << "ERROR a3 != a4" << std::endl;
    success &= (a5 == a6);
    if (! success) std::cerr << "ERROR a5 != a6" << std::endl;
    success &= (a7 == a8);
    if (! success) std::cerr << "ERROR a7 != a8" << std::endl;
    success &= (a9 == a10);
    if (! success) std::cerr << "ERROR a9 != a10" << std::endl;
    
    success &= (a1 == a3);
    if (! success) std::cerr << "ERROR a1 != a3" << std::endl;
    success &= (a5 == a7);
    if (! success) std::cerr << "ERROR a5 != a7" << std::endl;
    
    success &= (a1 == a5);
    if (! success) std::cerr << "ERROR a1 != a5" << std::endl;

    success &= (a1 == a9);
    if (! success) std::cerr << "ERROR a1 != a9" << std::endl;

#ifdef GIVARO_DEBUG
    if (! success)
        std::cerr << "Error: " << seed << std::endl;
    MemoryInfo.print(std::cerr) << std::endl;
#endif    

    ::Givaro::GivaroMain::End();
    return (! success);
}

