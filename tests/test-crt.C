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
#include <givaro/givmontg32.h>
#include <givaro/givzpz.h>
#include <givaro/givgfq.h>
#include <givaro/givcra.h>    // Chinese Remainder of two elements
#include <givaro/givrns.h>    // Chinese Remainder of an array of elements
#include <givaro/givrnsfixed.h>    // Chinese Remainder with fixed primes
#include <givaro/givtimer.h>
#include <givaro/givrandom.h>


typedef GFqDom<long> 		Field1;
typedef ZpzDom<Std16>           Field2;
typedef ZpzDom<Log16>           Field3;
typedef ZpzDom<Std32>  		Field4;
typedef ZpzDom<Std64>  		Field5;
typedef ZpzDom<Unsigned32>	Field6;
typedef Montgomery<Std32>       Field7;
typedef ZpzDom<Integer>         Field8;

template <typename Field>
Integer tmain(int argc, char ** argv, const GivRandom& generator) {
    typedef RNSsystem<Integer, Field >  CRTSystem;
    typedef typename CRTSystem::domains	Domains;
    typedef typename CRTSystem::array	Elements;
    typedef typename CRTSystem::ring	Ring;

    typedef RNSsystemFixed<Integer>  CRTSystemFixed;
    typedef typename RNSsystemFixed<Integer>::array  Prime_t;

    IntPrimeDom ID;
    Integer a( generator() >>(argc>2?atoi(argv[2]):17) ), M(1), b;

    Prime_t Primes( argc>1 ? atoi(argv[1]):15);
    Domains PrimeDoms( Primes.size() );
    Elements Moduli( Primes.size() );
    Prime_t ModuliInts( Primes.size() );

    typename Prime_t::iterator p = Primes.begin();
    typename Domains::iterator i = PrimeDoms.begin();
    typename Elements::iterator e = Moduli.begin();
    typename Prime_t::iterator m = ModuliInts.begin();
    for(; i != PrimeDoms.end(); ++i, ++e, ++p, ++m) {
        *i = Field( ID.nextprimein( a ) );
        *p = a;
//         i->random( generator, *e );
        i->init(*e,  generator() );
        *m = *e;
        M *= a;
    }

    CRTSystem CRT( PrimeDoms );

    Timer tim; tim.clear(); tim.start();
    CRT.RnsToRing( a, Moduli );
    tim.stop();
#ifdef GIVARO_DEBUG
    Field(PrimeDoms.front()).write( std::cerr << tim << " using ") << std::endl;
#endif

    if (PrimeDoms.size() < 50) {
        i = PrimeDoms.begin();
        e = Moduli.begin();
#ifdef GIVARO_DEBUG
        for( ; i != PrimeDoms.end(); ++i, ++e)
            i->write(std::cout << a << " mod " << i->characteristic() << " = ", *e) << ";" << std::endl;
#endif
    }

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

#ifdef GIVARO_DEBUG
    if (a != b)
        std::cerr << "incoherency between normal : " << a 
                  << " and fixed : " << b << std::endl;
#endif

   Elements Verifs( Primes.size() );
   CRT.RingToRns( Verifs, a );

#ifdef GIVARO_DEBUG
   typename Elements::const_iterator v = Verifs.begin();
   i = PrimeDoms.begin();
   e = Moduli.begin();
   for( ; i != PrimeDoms.end(); ++i, ++e, ++v)
       if (! i->areEqual(*e, *v) ) {
           i->write( std::cerr << "incoherency within ") << std::endl;
           break;
       }
#endif


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



int main(int argc, char ** argv) {
        // argv[1] : number of primes
        // argv[2] : 2^{32-j} is size of primes
        // argv[3] : seed for generator

    GivRandom seedor( argc>3 ? atoi(argv[3]):1234 );
    unsigned long seed = seedor.seed();


    Integer a1 = tmain<Field1>(argc, argv, GivRandom(seed));
    Integer a2 = tmain<Field2>(argc, argv, GivRandom(seed));
    Integer a3 = tmain<Field3>(argc, argv, GivRandom(seed));
    Integer a4 = tmain<Field4>(argc, argv, GivRandom(seed));
    Integer a5 = tmain<Field5>(argc, argv, GivRandom(seed));
    Integer a6 = tmain<Field6>(argc, argv, GivRandom(seed));
    Integer a7 = tmain<Field7>(argc, argv, GivRandom(seed));
    Integer a8 = tmain<Field8>(argc, argv, GivRandom(seed));

    bool success = true;
    success &= (a1 == a2); 
    if (! success) std::cerr << "ERROR a1 != a2" << std::endl;
    success &= (a3 == a4); 
    if (! success) std::cerr << "ERROR a3 != a4" << std::endl;
    success &= (a6 == a5); 
    if (! success) std::cerr << "ERROR a5 != a6" << std::endl;
    success &= (a7 == a8); 
    if (! success) std::cerr << "ERROR a7 != a8" << std::endl;
    success &= (a1 == a3); 
    if (! success) std::cerr << "ERROR a1 != a3" << std::endl;
    success &= (a5 == a7); 
    if (! success) std::cerr << "ERROR a5 != a7" << std::endl;
    success &= (a1 == a5); 
    if (! success) std::cerr << "ERROR a1 != a5" << std::endl;



#ifdef GIVARO_DEBUG
    if (! success)
        std::cerr << "Error: " << seed << std::endl;
#endif

    return (! success);
}
