// =================================================================== //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <16 Jun 15 16:05:18 Jean-Guillaume.Dumas@imag.fr>
// =================================================================== //

/** @file givintrsa.h
 * @ingroup integers
 * @brief RSA scheme.
 */

#ifndef __GIVARO_rsa_H
#define __GIVARO_rsa_H

#include <iostream>
#include "givaro/givinteger.h"
#include "givaro/givintprime.h"
#include "givaro/givintfactor.h"
#include "givaro/givrandom.h"


// k = 2^16 + 1, is prime
#define SIMPLE_EXPONENT (Element( (1<<16)+1 ))

// =================================================================== //
// RSA public-key cipher codes
// =================================================================== //

namespace Givaro {

    //! RSA domain
    template<class MyRandIter = GivRandom>
    class IntRSADom : public IntFactorDom<MyRandIter> {
    public:
        // JGD 19.02.2003 : Should work nicely, but produces:
        //                  "implicit typename is deprecated"
        // using IntFactorDom<MyRandIter>::Element;
        // using IntFactorDom<MyRandIter>::random_generator;
        typedef typename IntFactorDom<MyRandIter>::Element Element;
        typedef typename IntFactorDom<MyRandIter>::random_generator random_generator;

    private:
        Element _n, _e;
        Element _d;
        int64_t _lm;

    public:

        // =================================================================== //
        //! Constructors
        // =================================================================== //
        //@{
        IntRSADom(bool fi = false, MyRandIter g = MyRandIter() ) :
            IntFactorDom<MyRandIter>(g), _fast_impl(fi)
        {
            keys_gen(IntFactorDom<MyRandIter>::_g, 257, 255, _n, _e, _d);
            _lm = log(_n,1<<(8*sizeof(unsigned char)));
        }
        IntRSADom(const int64_t s, bool fi = false, MyRandIter g = MyRandIter() ) :
            IntFactorDom<MyRandIter>(g), _fast_impl(fi)
        {
            keys_gen(IntFactorDom<MyRandIter>::_g, (s>>1)-1, (s>>1)+1, _n, _e, _d);
            _lm = log(_n,1<<(8*sizeof(unsigned char)));
        }
        IntRSADom(const int64_t p, const int64_t q, bool fi = false, MyRandIter g = MyRandIter() ) :
            IntFactorDom<MyRandIter>(g), _fast_impl(fi)
        {
            keys_gen(IntFactorDom<MyRandIter>::_g, p, q, _n, _e, _d);
            _lm = log(_n,1<<(8*sizeof(unsigned char)));
        }
        IntRSADom(const Element& n, const Element& e, const Element& d) :
            _n(n), _e(e), _d(d), _lm(log(n,1<<(8*sizeof(unsigned char)))), _fast_impl( e == SIMPLE_EXPONENT )
        {}
        IntRSADom(const Element& n, const Element& e) :
            _n(n), _e(e), _d(0), _lm(log(n,1<<(8*sizeof(unsigned char)))), _fast_impl( e == SIMPLE_EXPONENT )
        {}
        //@}

        // =================================================================== //
        //! Accesses
        // =================================================================== //
        //@{
        const Element& getn() const { return _n; }
        const Element& gete() const { return _e; }
        const Element& getd() const { return _d; }
        //@}

        // =================================================================== //
        //! Text conversions
        // =================================================================== //
        //@{
        std::ostream& encipher(std::ostream&, std::istream&) const ;
        std::ostream& decipher(std::ostream&, std::istream&) ;
        //@}

        // =================================================================== //
        // Keys generation
        // public keys are n and e, the secret key is d.
        // ciphering is computing       : x^e mod m, in CBC mode
        // deciphering is computing     : b^d mod m
        // since for any x, x^(e.d) = x mod m
        // =================================================================== //


        // =================================================================== //
        /*! Strong Primes.
         * @bib
         * - J. Gordon, <i>Strong Primes Are Easy to Find</i>,  EUROCRYPT'84, LNCS 209.
         */
        // =================================================================== //
        Element& strong_prime(random_generator& g, int64_t psize, Element& p) const;

        // =================================================================== //
        /** Key gen.
         * Here m = p*q
         * p and q are prime numbers of respective sizes psize, qsize
         * Moreover p-1 and q-1 have one prime factor of respective size 2/3
         * since k.u = 1 mod (p-1)(q-1)
         */
        // =================================================================== //
        //@{
        void keys_gen(random_generator& g, int64_t psize, int64_t qsize, Element& n, Element& e, Element& d, Element& p, Element& q) const ;
        void keys_gen(random_generator& g, int64_t psize, int64_t qsize, Element& n, Element& e, Element& d) const ;
        //@}

        // =================================================================== //
        //! log[10]
        // =================================================================== //
        int64_t log(const Element& n, const int64_t = 10) const ;

        // =================================================================== //
        //! Text conversions
        // =================================================================== //
        std::ostream& ecriture_str(std::ostream&, const Element&) const ;
        std::ostream& ecriture_str_last(std::ostream&, const Element&) const ;
        std::ostream& ecriture_Int(std::ostream&, const Element&) const ;

        // =================================================================== //
        //! Breaking codes : finding u knowing only m an k ...
        // =================================================================== //
        Element& point_break(Element& u) ;

    protected:
        //! Fast implementation.
        //! Means simple enciphering key, and deciphering via chinese remaindering.
        //! @warning this means less security !
        bool _fast_impl;

    };

} // Givaro

#include "givaro/givintrsa.inl"

#endif // __GIVARO_rsa_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
