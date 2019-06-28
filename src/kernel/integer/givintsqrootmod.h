// ============================================================= //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <28 Jun 19 18:19:52 Jean-Guillaume.Dumas@imag.fr>
// Author : Yanis Linge adn Jean-Guillaume Dumas
// ============================================================= //


/** @file givintsqrootmod.h
 * @ingroup integers
 * @brief Modular square roots
 */

#ifndef __GIVARO_sqrtmod_H
#define __GIVARO_sqrtmod_H

#include <iostream>
#include "givaro/givinteger.h"
#include "givaro/givintprime.h"
#include "givaro/givintfactor.h"
#include "givaro/givintrns.h"
#include "givaro/givrandom.h"
#include "givaro/givpower.h"
#include <cmath>

namespace Givaro {

    //!  Modular square roots.
    template < class MyRandIter = GivRandom >
    class IntSqrtModDom:public IntFactorDom < MyRandIter > {
    public:
        typedef IntFactorDom < MyRandIter > Father_t;
        typedef typename IntFactorDom < MyRandIter >::Rep Rep;

        IntSqrtModDom (MyRandIter g = MyRandIter ()) : IntFactorDom < MyRandIter > (g) {}

        // ======================================================== //
        // Modular Square root functions
        // ======================================================== //
        Rep & sqrootmod (Rep & x, const Rep & a, const Rep & n) const ;

        // ======================================================== //
        // Brillhart decomposition of a prime
        // into a sum of squares
        // See [Note on representing a prime as a sum of two squares,
        //      J. Brillhart. Math. Comp. 26 (1972), 1011-1013 ]
        // ======================================================== //
        void Brillhart(Rep&, Rep&, const Rep&) const;

        // ======================================================== //
        // Element as a modular sum of squares
        // ======================================================== //
        void sumofsquaresmodprime(Rep&, Rep&, const Rep&, const Rep&) const;

        // ======================================================== //
        // Modular Square root sub-functions
        // ======================================================== //

        // p is supposed to be prime
        Rep & sqrootmodprime (Rep & x, const Rep & a, const Rep & p) const;

        // p is supposed to be prime, modulo is taken mod p^k
        Rep & sqrootmodprimepower (Rep & x, const Rep & a, const Rep & p, const uint64_t k, const Rep &pk) const;

        // modulo is taken mod 2^k
        Rep & sqrootmodpoweroftwo (Rep & x, const Rep & a,const uint64_t k, const Rep & pk) const;

    protected:

        // ======================================================== //
        // Linear update using only onemorelift
        // ======================================================== //

        // p is supposed to be prime and odd
        Rep & sqrootlinear (Rep & x, const Rep & a,const Rep & p,const uint64_t k) const;

        // result is modulo 2^{k+1}
        Rep & sqroottwolinear(Rep & x, const Rep & a,const uint64_t k) const;

        // ======================================================== //
        // Liftings
        // ======================================================== //

        // PRECONDITION: x^2 = a [p^k]
        // RETURNS: x s.t. x^2 = a [p^{2k}]
        Rep & sqroothensellift (Rep & x, const Rep & a, const Rep & p, const uint64_t k, const Rep & pk) const;

        // PRECONDITION: x^2 = a [p^k]
        // RETURNS: x s.t. x^2 = a [p^{k+1}]
        Rep & sqrootonemorelift (Rep & x, const Rep & a, const Rep & p, const uint64_t k, const Rep & pk) const;

        // PRECONDITION: x^2 = a [2^k], with k>=3, a and x are odd
        // RETURNS: x s.t. x^2 = a [2^{2k-2}]
        Rep & sqrootmodtwolift (Rep & x, const Rep & a, const uint64_t k, const Rep & pk) const;

    };

} // Givaro

#include "givaro/givintsqrootmod.inl"

#endif // __GIVARO_sqrtmod_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
