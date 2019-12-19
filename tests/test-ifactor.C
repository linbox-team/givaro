// Copyright(c)'1994-2010, 2011 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// written by BB
// see the COPYRIGHT file for more details.

//#include <iostream>
//#include <givaro/givinit.h>
#include <givaro/givintfactor.h>
//#include <givaro/givintprime.h>

#define NB_ITERS 15
//#define LOOPS 0

using namespace Givaro;

int test(const IntFactorDom<> & IP, const Integer & m)
{

    Integer f = 1 ;
    Integer r = 1 ;

    //IP.factor(f,m,LOOPS);
    IP.factor(f,m) ; // ne teste que Lenstra ou Pollard selon que que GIVARO_LENSTRA est définie ou non

    Integer::mod(r,m,f);
    if (r || ((f==1 || f==m) && !IP.local_prime(f)))
    {
#ifdef __GIVARO_DEBUG
        if (r)
            std::cout << "error : factor does not divide integer" << std::endl;
        else
            std::cout << "error : " << f << " is no good..." << std::endl;
#endif
        return -1 ;
    }


    IP.primefactor(f,m) ; // ne teste que Lenstra ou Pollard selon que que GIVARO_LENSTRA est définie ou non
    if (r || !IP.local_prime(f))
    {
#ifdef __GIVARO_DEBUG
        if (r)
            std::cout << "error : factor does not divide integer" << std::endl;
        else
            std::cout << "error : " << f << " is no good..." << std::endl;
#endif
        return -1 ;
    }

    return 0 ;
}

int main()
{
    IntFactorDom<> IP;
    Integer m;
    m.seeding();
    int err = 0 ;

    long int a = 3;
    long int b = 50 ;
    for (size_t i = 0 ; i < NB_ITERS ; ++i)
    {
        //if (!(i%25)) std::cout << i << "..." ;
        m = Integer::random_between(a,b);
        err = test(IP,m);
        if (err) break ;
    }
    if (err) return err ;

    a = 1;
    b = 30 ;
    for (size_t i = 0 ; i < NB_ITERS ; ++i)
    {
        m = Integer::random_between(a,b);
        err = test(IP,m);
        if (err) break ;
    }
    if (err) return err ;

    // harder :
    Integer p,q ;
    a = 19 ;
    b = 20 ;
    for (size_t i = 0 ; i < NB_ITERS/2 ; ++i)
    {
        p = Integer::random_between(a,b);
        IP.nextprimein(p);
        q = Integer::random_between(a,b);
        IP.nextprimein(q);

        m = p*q ;
        err = test(IP,m);
        if (err) break ;
    }
    if (err) return err ;

    a = 25 ;
    b = 26 ;
    for (size_t i = 0 ; i < NB_ITERS/2 ; ++i)
    {
        p = Integer::random_between(a,b);
        IP.nextprimein(p);
        q = Integer::random_between(a,b);
        IP.nextprimein(q);

        m = p*q ;
        err = test(IP,m);
        if (err) break ;
    }
    if (err) return err ;

#ifdef __GIVARO_DEBUG
    std::cout << "success" <<std::endl;
#endif
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
