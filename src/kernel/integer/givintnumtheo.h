// =================================================================== //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
//
// Time-stamp: <28 Jun 19 16:49:39 Jean-Guillaume.Dumas@imag.fr>
// =================================================================== //

/*! @file givintnumtheo.h
 * @ingroup integers
 * @brief num theory.
 * - Euler's phi function
 * -  Primitive roots
 * - RSA scheme.
 * .
 */

#ifndef __GIVARO_numtheory_H
#define __GIVARO_numtheory_H

#include <iostream>
#include "givaro/givinteger.h"
#include "givaro/givintprime.h"
#include "givaro/givintfactor.h"
#include "givaro/givrandom.h"
#include "givaro/udl.h"

namespace Givaro {
    //! Num theory Domain.
    template<class MyRandIter = GivRandom>
    class IntNumTheoDom : public IntFactorDom<MyRandIter> {
    public:
        typedef IntFactorDom<MyRandIter> Father_t;
        typedef typename IntFactorDom<MyRandIter>::Rep Rep;
        IntNumTheoDom(MyRandIter g = MyRandIter())
        :  IntFactorDom<MyRandIter>(g) {}

        // =================================================================== //
        //! Euler's phi function
        // =================================================================== //
        template <template <class, class> class Container, template<class> class Alloc>
        Rep& phi(Rep& res, const Container<Rep, Alloc<Rep> >& Lf, const Rep& n) const ;

        Rep& phi(Rep& r, const Rep& n) const ;

        // =================================================================== //
        //! Primitive Root
        // =================================================================== //
        Rep& prim_root(Rep&, const Rep&) const ;
        Rep& prim_root(Rep&, uint64_t&, const Rep&) const ;
        Rep& prim_root_of_prime(Rep&, const Rep&) const ;
        template<class Array> Rep& prim_root_of_prime(Rep& A, const Array& Lf, const Rep& phin, const Rep& n) const ;

        /**  Polynomial-time generation of primitive roots.
         *  L is number of loops of Pollard partial factorization of n-1
         *  10,000,000 gives at least 1-2^{-40} probability of success
         *  [Dubrois & Dumas, Industrial-strength primitive roots]
         *  Returns the probable primitive root and the probability of error.
         */
        Rep& probable_prim_root(Rep&, double&, const Rep& n, const uint64_t L = 10000000_ui64) const;

        //!  Here L is computed so that the error is close to epsilon
        Rep& probable_prim_root(Rep&, double&, const Rep& n, const double epsilon) const;

        Rep& lowest_prim_root(Rep&, const Rep&) const ;
        bool is_prim_root(const Rep&, const Rep&) const ;
        Rep& order(Rep&, const Rep&, const Rep&) const ;
        bool isorder(const Rep&, const Rep&, const Rep&) const ;

        // =================================================================== //
        /** Generalization of primitive roots for any modulus
         * Primitive means maximal order
         *    Primitive Element, Primitive invertible
         *    Both functions coincide except for m=8
         *
         * Lambda Function : maximal orbit size
         *    lambda : Order of a primitive Element
         *    lambda_inv : Order of an invertible primitive Element
         *    Both functions coincide except for m=8
         */
        // =================================================================== //
        Rep& prim_inv(Rep & , const Rep&) const ;
        Rep& prim_elem(Rep & , const Rep&) const ;
    private:
        Rep& prim_base(Rep & , const Rep&) const ;
        Rep& lambda_base(Rep & , const Rep&) const ;
    public:
        Rep& lambda_primpow(Rep & , const Rep&, uint64_t) const ;
        Rep& lambda_inv_primpow(Rep & , const Rep&, uint64_t) const ;
        Rep& lambda(Rep & , const Rep&) const ;
        Rep& lambda_inv(Rep & , const Rep&) const ;

        // =================================================================== //
        //! Möbius function
        // =================================================================== //
        template< template<class, class> class Container, template <class> class Alloc>
        short mobius(const Container<Rep, Alloc<Rep> >& lpow) const ;

        //! Möbius function
        short mobius(const Rep& a) const;
    };

} // givaro
#include "givaro/givintnumtheo.inl"

#endif // __GIVARO_numtheory_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
