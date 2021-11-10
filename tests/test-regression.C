// Copyright(c)'1994-2025 by The Givaro group
// This file is part of Givaro.
// BB: adapted,enhanced from examples/FiniteFields/ff-arith.C
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

#include <iostream>

#ifdef __GIVARO_DEBUG
#define GIVARO_RATRECON_DEBUG
#endif

#include <givaro/givinit.h>
#include <givaro/qfield.h>
using namespace Givaro;


bool testRationalDenom() {

    bool pass = true;

    QField<Rational> QQ;

    Rational a(-1,1), b(1,-1), c;
    std::clog << "-1: " << a << std::endl;
    std::clog << "-1: " << b << std::endl;
    pass = pass && (a == QQ.mOne);
    pass = pass && (b == QQ.mOne);

    QQ.invin(a);
    std::clog << "-1: " << a << std::endl;
    pass = pass && (a == QQ.mOne);

    QQ.inv(a, b);
    std::clog << "-1: " << a << std::endl;
    pass = pass && (a == QQ.mOne);

    QQ.negin(a);
    std::clog << "1: " << a << std::endl;
    pass = pass && (a == QQ.one);

    QQ.neg(a,b);
    std::clog << "1: " << a << std::endl;
    pass = pass && (a == QQ.one);

    QQ.add(c,a,b);
    std::clog << "0: " << c << std::endl;
    pass = pass && isZero(c);

    std::clog << "[TRD] " << (pass?"PASSED.":"FAILED.") << std::endl;

    return pass;
}


#include <givaro/modular.h>
bool testFieldInit() {
    size_t p = 19;
    Givaro::Modular<int64_t> F1(p);

        // This is a compilation test,
        // if it compiles, everything's fine
    F1.write(std::clog << "[TFI] PASSED: ") << std::endl;
    return true;
}

#include <givaro/givintprime.h>
bool testPrevNextPrime() {

    bool pass = true;

    IntPrimeDom IP;

    Integer a(1),b(1),p,q,r;
    a<<=125 ;
    b<<=136 ;

    p = Integer::random_between(a,b);
    std::clog << p << std::endl;

    IP.nextprimein(p);
    IP.nextprimein(p);
    IP.prevprime(q,p);
    IP.prevprimein(q);
    IP.nextprimein(q);
    IP.nextprime(r,q);
    pass = IP.areEqual(p,r);

    std::clog << "[PRI] " << (pass?"PASSED.":"FAILED.") << std::endl;

    return pass;
}


int main(int argc, char ** argv)
{
    bool pass = true;

    const int seed = int (argc>1?atoi(argv[1]):BaseTimer::seed());
#ifdef __GIVARO_DEBUG
    GivaroMain::DisplayVersion(std::clog);
    std::clog << "seed: " << seed << std::endl;
#endif
    Integer::seeding((uint64_t)seed);


    pass = pass && testRationalDenom  ();
    pass = pass && testFieldInit  ();
    pass = pass && testPrevNextPrime  ();

    return pass ? 0 : -1;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
