// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
/*! @file examples/Polynomial/PolynomialCRT.C
 * @ingroup examples
 * @ingroup polynomials
 * @example examples/Polynomial/PolynomialCRT.C
 * @brief NO DOC
 */

#include <iostream>
#include <algorithm>
#include <givaro/givtimer.h>
#include <givaro/givpoly1crt.h>
#include <givaro/givintprime.h>
#include <givaro/montgomery.h>
#include <givaro/extension.h>
#include <givaro/modular.h>
#include <givaro/gfq.h>
#include <givaro/chineseremainder.h>    // Chinese Remainder of two elements
#include <givaro/givrns.h>    // Chinese Remainder of an array of elements
#include <givaro/givrandom.h>
#include <givaro/qfield.h>


using namespace Givaro;


typedef GFqDom<int64_t> 	Field1;
typedef Modular<int16_t>    Field2;
typedef Modular<Log16>      Field3;
typedef Modular<int32_t>  	Field4;
typedef Modular<int64_t>  	Field5;
typedef Modular<uint32_t>	Field6;
typedef Montgomery<int32_t> Field7;
typedef QField<Rational> 	Field8;
typedef Modular<Integer>    Field9;


typedef Extension<> 		Field10;

template <typename Field>
bool tmain(int argc, char ** argv, GivRandom& generator) {
    bool pass = true;
    typedef Poly1CRT< Field >  CRTSystem;
    typedef typename CRTSystem::Element	Poly;
    //     typedef typename CRTSystem::Type_t	Scal;
    //     typedef typename CRTSystem::array_E	VPoly;
    typedef typename CRTSystem::array_T	VScal;

    IntPrimeDom ID;
    Integer a( generator() >>(argc>2?atoi(argv[2]):17) );
    Field F(ID.nextprimein( a ));

    VScal Primes( argc>1 ? (size_t)atoi(argv[1]):15);
    VScal Moduli( Primes.size() );

    typename VScal::iterator i = Primes.begin();
    typename VScal::iterator e = Moduli.begin();
    for(; i != Primes.end(); ++i, ++e) {
        do {
            F.init(*i, generator());
        } while ( (std::find(Primes.begin(), i, *i) != i) || (F.isZero(*i))) ;

        F.init(*e, generator());
    }
    //    for(typename VScal::const_iterator it=Primes.begin(); it!=Primes.end();++it)
    //        F.write(std::cout, *it) << std::endl;
    //    for(typename VScal::const_iterator it=Moduli.begin(); it!=Moduli.end();++it)
    //        F.write(std::cout, *it) << std::endl;


    CRTSystem CRT( F, Primes, "X" );
    Poly res;

    Timer tim; tim.clear(); tim.start();
    CRT.RnsToRing( res, Moduli );
    tim.stop();
    F.write( std::cerr << tim << " using ") << std::endl;

    if (Primes.size() < 14) {
        i = Primes.begin();
        e = Moduli.begin();
        for( ; i != Primes.end(); ++i, ++e)
            if (F.characteristic()>0)
                F.write(CRT.getpolydom().write(F.write(std::cout << "subs(X=", *i) << ",", res) << ") mod " << F.characteristic() << " = ", *e) << ';' << std::endl;
            else
                F.write(CRT.getpolydom().write(F.write(std::cout << "subs(X=", *i) << ",", res) << ") = ", *e) << ';' << std::endl;
    }


    VScal Verifs( Primes.size() );
    CRT.RingToRns( Verifs, res );
    typename VScal::const_iterator v = Verifs.begin();
    e = Moduli.begin();
    for( ; e != Moduli.end(); ++e, ++v)
        if (! F.areEqual(*e, *v) ) {
            F.write(std::cerr << "incoherency within ") << std::endl;
            F.write(std::cerr << "e: ", *e ) << std::endl;
            F.write(std::cerr << "v: ", *v ) << std::endl;
            pass = false;
            break;
        }

    CRT.getpolydom().random(generator, res, Degree((int64_t)Primes.size()-1));
    CRT.RingToRns( Verifs, res );
    Poly nres;

    tim.clear(); tim.start();
    CRT.RnsToRing( nres, Verifs );
    tim.stop();
    if (! CRT.getpolydom().areEqual(res,nres) ) {
        CRT.getpolydom().write(std::cerr << "incoherency within ") << std::endl;
        CRT.getpolydom().write(std::cerr << "r: ", res ) << std::endl;
        CRT.getpolydom().write(std::cerr << "n: ", nres ) << std::endl;
        pass = false;
    }
    F.write( std::cerr << tim << " using ") << std::endl;


    return pass;
}
template <typename Field>
bool tmainext(int argc, char ** argv, GivRandom& generator) {
    bool pass = true;
    typedef Poly1CRT< Field >  CRTSystem;
    typedef typename CRTSystem::Element	Poly;
    //     typedef typename CRTSystem::Type_t	Scal;
    //     typedef typename CRTSystem::array_E	VPoly;
    typedef typename CRTSystem::array_T	VScal;

    IntPrimeDom ID;
    Integer a( generator() >>(argc>2?atoi(argv[2]):17) );
    Field F(ID.nextprimein( a ),2);

    VScal Primes( argc>1 ? (size_t)atoi(argv[1]):15);
    VScal Moduli( Primes.size() );

    typename VScal::iterator i = Primes.begin();
    typename VScal::iterator e = Moduli.begin();
    for(; i != Primes.end(); ++i, ++e) {
        F.init(*i); F.init(*e);
        do {
            F.random(generator,*i);
        } while ( (std::find(Primes.begin(), i, *i) != i) || (F.isZero(*i))) ;

        F.random(generator,*e);
    }
    for(typename VScal::const_iterator it=Primes.begin(); it!=Primes.end();++it)
        F.write(std::cout, *it) << std::endl;
    for(typename VScal::const_iterator it=Moduli.begin(); it!=Moduli.end();++it)
        F.write(std::cout, *it) << std::endl;


    CRTSystem CRT( F, Primes, "X" );
    Poly res;

    Timer tim; tim.clear(); tim.start();
    CRT.RnsToRing( res, Moduli );
    tim.stop();
    F.write( std::cerr << tim << " using ") << std::endl;

    if (Primes.size() < 14) {
        i = Primes.begin();
        e = Moduli.begin();
        for( ; i != Primes.end(); ++i, ++e)
            if (F.characteristic()>0)
                F.write(CRT.getpolydom().write(F.write(std::cout << "subs(X=", *i) << ",", res) << ") mod " << F.characteristic() << " = ", *e) << ';' << std::endl;
            else
                F.write(CRT.getpolydom().write(F.write(std::cout << "subs(X=", *i) << ",", res) << ") = ", *e) << ';' << std::endl;
    }


    VScal Verifs( Primes.size() );
    CRT.RingToRns( Verifs, res );
    typename VScal::const_iterator v = Verifs.begin();
    e = Moduli.begin();
    for( ; e != Moduli.end(); ++e, ++v)
        if (! F.areEqual(*e, *v) ) {
            F.write(std::cerr << "incoherency within ") << std::endl;
            F.write(std::cerr << "e: ", *e ) << std::endl;
            F.write(std::cerr << "v: ", *v ) << std::endl;
            pass = false;
            break;
        }

    CRT.getpolydom().random(generator, res, Degree((int64_t)Primes.size()-1));
    CRT.RingToRns( Verifs, res );
    Poly nres;

    tim.clear(); tim.start();
    CRT.RnsToRing( nres, Verifs );
    tim.stop();
    if (! CRT.getpolydom().areEqual(res,nres) ) {
        CRT.getpolydom().write(std::cerr << "incoherency within ") << std::endl;
        CRT.getpolydom().write(std::cerr << "r: ", res ) << std::endl;
        CRT.getpolydom().write(std::cerr << "n: ", nres ) << std::endl;
        pass = false;
    }
    F.write( std::cerr << tim << " using ") << std::endl;


    return pass;
}



int main(int argc, char ** argv) {
    // argv[1] : number of primes
    // argv[2] : 2^{32-j} is size of primes
    // argv[3] : seed for generator

    GivRandom seedor( argc>3 ? (uint64_t)atoi(argv[3]): (uint64_t)BaseTimer::seed() );
    uint64_t seed = seedor.seed();
    std::cerr << "seed: " << seed << std::endl;

    Integer::seeding(seed);

    return
    tmain<Field1>(argc, argv, *( new GivRandom(seed)))
    && tmain<Field2>(argc, argv, *( new GivRandom(seed)))
    && tmain<Field3>(argc, argv, *( new GivRandom(seed)))
    && tmain<Field4>(argc, argv, *( new GivRandom(seed)))
    && tmain<Field5>(argc, argv, *( new GivRandom(seed)))
    && tmain<Field6>(argc, argv, *( new GivRandom(seed)))
    && tmain<Field7>(argc, argv, *( new GivRandom(seed)))
    && tmain<Field8>(argc, argv, *( new GivRandom(seed)))
    && tmain<Field9>(argc, argv, *( new GivRandom(seed)))
    && tmainext<Field10>(argc, argv, *( new GivRandom(seed)))
    ;

}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
