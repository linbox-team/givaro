// ============================================================= //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <19 Jan 11 15:53:22 Jean-Guillaume.Dumas@imag.fr>
// Author : Yanis Linge
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
#include <math.h>

namespace Givaro {

	//!  Modular square roots.
template < class RandIter = GivRandom >
class IntSqrtModDom:public IntFactorDom < RandIter > {
public:
    typedef IntFactorDom < RandIter > Father_t;
    typedef typename IntFactorDom < RandIter >::Rep Rep;

    IntSqrtModDom (RandIter g = RandIter ()) : IntFactorDom < RandIter > (g) {}

        // ======================================================== //
        // Modular Square root functions
        // ======================================================== //
    Rep & sqrootmod (Rep & x, const Rep & a, const Rep & n) const {
        std::vector < Rep > Lf;
        std::vector < unsigned long > Le;
        Father_t::set (Lf, Le, n);

        typename std::vector < Rep >::const_iterator Lf_iter = Lf.begin ();
        typename std::vector < unsigned long >::const_iterator Le_iter = Le.begin ();

        std::vector < Rep > roots;
        Rep tmp;

            // Build prime powers
        std::vector < Rep > Pe (Lf.size ());
        typename std::vector < Rep >::iterator Pe_iter = Pe.begin ();
        for (; Pe_iter != Pe.end (); ++Pe_iter, ++Lf_iter, ++Le_iter)
			dom_power (*Pe_iter, *Lf_iter, (long)*Le_iter, *this);

        Lf_iter = Lf.begin ();
        Le_iter = Le.begin ();
        Pe_iter = Pe.begin ();

            // roots mod powers of primes
        for (; Lf_iter != Lf.end (); ++Lf_iter, ++Le_iter, ++Pe_iter){
                // root mod a power of 2
            if (*Lf_iter == 2UL){
                roots.push_back (
                    this->sqrootmodpoweroftwo (tmp, a, *Le_iter, *Pe_iter));
            } else {
                roots.push_back (
                    this->sqrootmodprimepower (tmp, a, *Lf_iter, *Le_iter, *Pe_iter));
            }
        }

            // Chinese Remaindering
        IntRNSsystem < std::vector, std::allocator > RNs (Pe);

        RNs.RnsToRing (x, roots);
        x = (x<0?-x:x);
        return x;
    }

        // ======================================================== //
        // Modular Square root sub-functions
        // ======================================================== //

        // p is supposed to be prime
    Rep & sqrootmodprime (Rep & x, const Rep & a, const Rep & p) const;

        // p is supposed to be prime, modulo is taken mod p^k
    Rep & sqrootmodprimepower (Rep & x, const Rep & a, const Rep & p, const unsigned long k, const Rep &pk) const;

        // modulo is taken mod 2^k
    Rep & sqrootmodpoweroftwo (Rep & x, const Rep & a,const unsigned long k, const Rep & pk) const;

protected:

        // ======================================================== //
        // Linear update using only onemorelift
        // ======================================================== //

        // p is supposed to be prime and odd
    Rep & sqrootlinear (Rep & x, const Rep & a,const Rep & p,const unsigned long k) const;

        // result is modulo 2^{k+1}
    Rep & sqroottwolinear(Rep & x, const Rep & a,const unsigned long k) const;

        // ======================================================== //
        // Liftings
        // ======================================================== //

        // PRECONDITION: x^2 = a [p^k]
        // RETURNS: x s.t. x^2 = a [p^{2k}]
    Rep & sqroothensellift (Rep & x, const Rep & a, const Rep & p, const unsigned long k, const Rep & pk) const;

        // PRECONDITION: x^2 = a [p^k]
        // RETURNS: x s.t. x^2 = a [p^{k+1}]
    Rep & sqrootonemorelift (Rep & x, const Rep & a, const Rep & p, const unsigned long k, const Rep & pk) const;

        // PRECONDITION: x^2 = a [2^k], with k>=3, a and x are odd
        // RETURNS: x s.t. x^2 = a [2^{2k-2}]
    Rep & sqrootmodtwolift (Rep & x, const Rep & a, const unsigned long k, const Rep & pk) const;

};

} // Givaro

#include "givaro/givintsqrootmod.inl"

#endif // __GIVARO_sqrtmod_H
