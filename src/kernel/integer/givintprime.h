// =================================================================== //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <10 Nov 21 16:20:13 Jean-Guillaume.Dumas@imag.fr>
// =================================================================== //


/*! @file givintprime.h
 * @ingroup integers
 * @brief primes
 * - Prime numbers
 * - Modular powering,
 * - Fermat numbers,
 * - Primality tests
 * - Factorization : (There are parameters to fix)
 * .
 */
#ifndef __GIVARO_integers_prime_H
#define __GIVARO_integers_prime_H

#include "givaro/givinteger.h"

namespace Givaro {

    // =================================================================== //
    //! Fermat numbers
    // =================================================================== //
    class FermatDom : public IntegerDom {
    public:
        FermatDom() : IntegerDom() {}
        Rep& fermat (Rep&, const size_t)  const ;
        bool pepin (const size_t) const ;
    private:
        bool pepin (const Integer&) const ;
    };


    // =================================================================== //
    // Primality tests and factorization algorithms
    // =================================================================== //

    // Those macros are parameters to fix

    // primes known
    // first array
#define LOGMAX 3512
#define TABMAX 32768
    // second array
#define LOGMAX2 3031
#define TABMAX2 65536
    // Bounds between big and small
#define BOUNDARY_isprime TABMAX
#define BOUNDARY_2_isprime TABMAX2

#define GIVARO_ISLT(a,b) ((a)<(b))
#define GIVARO_ISLEQ(a,b) ((a)<=(b))
#define GIVARO_ISGT(a,b) ((a)>(b))
#define GIVARO_ISGEQ(a,b) ((a)>=(b))


    // =================================================================== //
    //! Primality tests
    // =================================================================== //
    class IntPrimeDom : public IntegerDom {
    public:
        IntPrimeDom() :  IntegerDom() {}

        int isprime(const Rep& n, int r=_GIVARO_ISPRIMETESTS_) const
        {
            /*
               return probab_prime(n);
               */
            //             return ((n)<BOUNDARY_isprime ?  isprime_Tabule(n) :
            //                     (n)<BOUNDARY_2_isprime ? isprime_Tabule2(n) :
            //                     probab_prime(n));
            int64_t l;
            return int32_t (int32_t(GIVARO_ISLT(n,BOUNDARY_isprime) ?  isprime_Tabule((int32_t)convert(l,n)):
                                    GIVARO_ISLT(n,BOUNDARY_2_isprime) ? isprime_Tabule2((int32_t)convert(l,n)):
                                    local_prime(n,r)));
        }

        // if p is a prime power, p = r^return
        // else return is 0 and r is undefined
        unsigned int isprimepower(Rep&, const Rep&) const ;

        template<class MyRandIter>
        unsigned int Miller(MyRandIter& g, const Rep& n=_GIVARO_ISPRIMETESTS_) const  ;

        template<class MyRandIter>
        Rep& test_Lehmann(MyRandIter& g, Rep&, const Rep& n=_GIVARO_ISPRIMETESTS_) const  ;

        template<class MyRandIter>
        int Lehmann(MyRandIter& g, const Rep& n=_GIVARO_ISPRIMETESTS_)  const ;

        int isprime_Tabule(const int n) const ;
        int isprime_Tabule2(const int n) const ;

        Rep& nextprime(Rep&, const Rep&, int r=_GIVARO_ISPRIMETESTS_) const ;
        Rep& prevprime(Rep&, const Rep&, int r=_GIVARO_ISPRIMETESTS_) const ;
        Rep& nextprimein(Rep&, int r=_GIVARO_ISPRIMETESTS_) const ;
        Rep& prevprimein(Rep&, int r=_GIVARO_ISPRIMETESTS_) const ;


        // Using Integer
        int local_prime(const Rep& n, int r=_GIVARO_ISPRIMETESTS_) const
        {
            return Protected::probab_prime(n,r);
        }


    private:
        static int IP[LOGMAX+5];  // -- table for Tabule
        static const int * TP;    // -- shifted table
        static int IP2[LOGMAX2+5]; // -- table for Tabule2
        static const int * TP2;    // -- shifted table
#if 0
        static int Tabule2(const Integer& p) ;
        static int Tabule(const Integer& p) ;
        static int _memTab2[LOGMAX2+5];   // -- table for Tabule2
        static const int* _Tab2; // -- shifted _memTabule2
        static int _memTab[];    // -- table for Tabule
        static const int* _Tab;  // -- shifted _memTabule
#endif
    };

} // Givaro
#include "givaro/givintprime.inl"
#endif // __GIVARO_integers_prime_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
