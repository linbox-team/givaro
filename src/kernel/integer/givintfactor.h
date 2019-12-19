// =================================================================== //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Needs Container structures : stl ones for instance
// Time-stamp: <16 Jun 15 16:05:40 Jean-Guillaume.Dumas@imag.fr>
// =================================================================== //

/*! @file givintfactor.h
 * @ingroup integers
 * @brief factorisation
 *- Prime numbers
 * - Factor sets :
 * - Pollard's rho method for factorization
 * - Elliptic curves factorization by Lenstra
 * .
 */

#ifndef __GIVARO_factorisation_H
#define __GIVARO_factorisation_H

#include <iostream>
#include "givaro/givinteger.h"
#include "givaro/givintprime.h"
#include "givaro/givrandom.h"

// #define BOUNDARY_factor TABMAX2

#define factor_first_primes(tmp,n) (tmp = isZero(mod(tmp,n,23))?23:( isZero(mod(tmp,n,19))?19:( isZero(mod(tmp,n,17))?17:  (isZero(mod(tmp,n,2))?2:( isZero(mod(tmp,n,3))?3:( isZero(mod(tmp,n,5))?5:( isZero(mod(tmp,n,7))?7: ( isZero(mod(tmp,n,11))?11:13 ))))))))

#define factor_second_primes(tmp,n) (tmp = isZero(mod(tmp,n,31))?31:( isZero(mod(tmp,n,29))?29: ( isZero(mod(tmp,n,37))?37: ( isZero(mod(tmp,n,41))?41:( isZero(mod(tmp,n,43))?43:  ( isZero(mod(tmp,n,71))?71:( isZero(mod(tmp,n,67))?67:( isZero(mod(tmp,n,61))?61:( isZero(mod(tmp,n,59))?59: ( isZero(mod(tmp,n,53))?53:( isZero(mod(tmp,n,47))?47: ( isZero(mod(tmp,n,97))?97: ( isZero(mod(tmp,n,89))?89:( isZero(mod(tmp,n,83))?83:( isZero(mod(tmp,n,79))?79:73)))))))))))))))

namespace Givaro {

    // =================================================================== //
    // Set or Container of divisors, factors.
    // =================================================================== //

    //! Integer Factor Domain.
    template<class MyRandIter = GivRandom>
    class IntFactorDom : public IntPrimeDom {
    private:
        // 2*3*5*7*11*13*17*19*23
        const int PROD_first_primes;
        // 29*31*37*41*43*47*53*59*61*67*71*73*79*83*89*97
        const Rep PROD_second_primes;
    protected:
        MyRandIter _g;

    public:
        typedef MyRandIter random_generator;

        IntFactorDom(MyRandIter g = MyRandIter()) :
            IntPrimeDom(),
            PROD_first_primes(223092870),
            PROD_second_primes("10334565887047481278774629361"),
            _g(g)
        {
        }

        //  loops defaulted to 0 forces Pollard's factorization to
        //  be complete
        Rep& factor(Rep& r, const Rep& n, unsigned long loops = 0) const
        {
            if (isOne(gcd(r,n,PROD_first_primes)))
                if (isOne(gcd(r,n,PROD_second_primes))) {
#ifdef GIVARO_LENSTRA
                    return Lenstra((const MyRandIter&)_g, r, n);
#else
                    return Pollard((const MyRandIter&)_g, r, n, loops);
#endif
                } else
                    return factor_second_primes(r,n);
                else
                    return factor_first_primes(r,n);
        }

        //  Factors are checked for primality
        Rep& iffactorprime (Rep& r, const Rep& n, unsigned long loops = 0) const
        {
            if (factor(r, n, loops) != 1) {
                if (! isprime(r,_GIVARO_ISPRIMETESTS_) ) {
                    Rep nn = r; factor(r,nn, loops);
                }
                while (! isprime(r,_GIVARO_ISPRIMETESTS_) ) {
                    Rep nn = r;
                    if (isOne(gcd(r,nn,PROD_first_primes))) {
                        if (isOne(gcd(r,nn,PROD_second_primes))) {
                            Pollard((const MyRandIter&)_g, r, nn, loops);
                        } else {
                            factor_second_primes(r,nn);
                        }
                    } else {
                        factor_first_primes(r,nn);
                    }
                    if (r == nn) {
                        Lenstra((const MyRandIter&)_g, r, nn) ;
                        break; // In case Lenstra fails also
                    }
                }
            }
            return r;
        }

        Rep& primefactor(Rep& r, const Rep& n) const {
            while ((iffactorprime(r,n,0) == 1) && (! isprime(n, _GIVARO_ISPRIMETESTS_)) ) {}
            return r;
        }


        /// Factors with primes
        //  loops defaulted to 0 forces factorization to be complete
        //  otherwise returns if factorization is complete or not
        //  Factors are checked for primality
        template<class Container1, class Container2> bool set
        ( Container1& setint, Container2& setpwd,  const Rep& a, unsigned long loops = 0) const ;
        ///
        template<class Container> void set( Container&,  const Rep&) const ;
        ///
        template<class Container> void Erathostene(Container&, const Rep&) const ;
        ///
        template<class Container, class Cont2, class Cont3> Container& divisors(Container& L, const Cont2& Lf, const Cont3& Le)  const;
        template<class Container> Container& divisors(Container&, const Rep& ) const ;

        /// returns a small factor
        Rep& Erathostene(Rep&,  const Rep& p ) const ;

        // Pollard with a bound on the number of loops
        // Bound 0 is without bound
        Rep& Pollard(const MyRandIter&, Rep&, const Rep& n, unsigned long threshold = 0) const ;
        // returns a factor by Lenstra's elliptic curves method
        Rep& Lenstra(const MyRandIter&, Rep&, const Rep& n, const Rep& B1 = 10000000, const unsigned long curves = 30) const ;

        std::ostream& write(std::ostream& o, const Rep& n) const;
        template<class Array> std::ostream& write(std::ostream& o, Array&, const Rep& n) const;


    private:
        // Those are parameters for Pollard's algorithms
        // Pollard_fctin : must be somewhat a "random" function in Z/nZ
        // Pollard_cst can be a slight alternative for the Pfct x^2+1
#ifndef Pollard_cst
#define Pollard_cst 1
#endif

        using IntPrimeDom::write; // to silence off clang warning about hiding overloaded function

        Rep& Pollard_fctin(Rep & x, const Rep& n) const
        {
            mulin(x,x);
            addin(x,Pollard_cst);
            return modin(x,n);
        }

    };

} // Givaro

#include "givaro/givintfactor.inl"

#endif // __GIVARO_factorisation_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
