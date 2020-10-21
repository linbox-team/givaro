// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: JG Dumas
//          A. Breust (adapted)
// ========================================================================

/*! @file givzpzInt.h
 * @ingroup zpz
 *  @brief Arithmetic on Z/pZ, with p a prime number in arbitrary precision.
 */

#ifndef __GIVARO_zpz_gen_H
#define __GIVARO_zpz_gen_H

#include "givaro/givbasictype.h"
#include "givaro/giverror.h"
#include "givaro/givinteger.h"
#include "givaro/givcaster.h"
#include "givaro/givranditer.h"
#include "givaro/givtypestring.h"
#include "givaro/ring-interface.h"
#include "givaro/modular-implem.h"
#include "givaro/modular-general.h"

namespace Givaro
{

    /*! @brief This class implement the standard arithmetic with Modulo Elements.
     * - The representation of an integer a in Zpz is the value a % p
     * .
     */
    template<typename IntType, typename _Compute_t, typename Enable>
    class Modular : public Modular_implem<IntType, _Compute_t, IntType>
    {
    public:

        // ----- Exported Types and constantes
        using Storage_t = IntType;
        using Residu_t = IntType;
        using Compute_t = _Compute_t;

        using Self_t = Modular<Storage_t, Compute_t>;
        using Parent_t = Modular_implem<Storage_t, Compute_t, Residu_t>;

        using Element = Storage_t;

        // ----- Constructors

        using Modular_implem<Storage_t, Compute_t, Residu_t>::Modular_implem;
        using Parent_t::_p;
        using Parent_t::_pc;
        using Parent_t::zero;
        using Parent_t::one;
        using Parent_t::mOne;

        // ----- Initialisation
        Element& init (Element& x) const
        { return x = 0; }
        Element& init (Element& x, const Integer& y) const final
        { x = y % _p; return reduce(x); }
        template<typename T> Element& init(Element& r, const T& a) const
        { r = Caster<Element>(a); return reduce(r); }


        Element& reduce (Element& x, const Element& y) const
        { x = y % _p; if (x < 0) x += _p; return x; }
        Element& reduce (Element& x) const
        { x %= _p; if (x < 0) x += _p; return x; }

        // ----- Classic arithmetic
        Element& mul(Element& r, const Element& a, const Element& b) const;
        Element& div(Element& r, const Element& a, const Element& b) const;
        Element& add(Element& r, const Element& a, const Element& b) const;
        Element& sub(Element& r, const Element& a, const Element& b) const;
        Element& neg(Element& r, const Element& a) const;
        Element& inv(Element& r, const Element& a) const;

        Element& mulin(Element& r, const Element& a) const;
        Element& divin(Element& r, const Element& a) const;
        Element& addin(Element& r, const Element& a) const;
        Element& subin(Element& r, const Element& a) const;
        Element& negin(Element& r) const;
        Element& invin(Element& r) const;

        // -- axpy:   r <- a * x + y
        // -- axpyin: r <- a * x + r
        Element& axpy  (Element& r, const Element& a, const Element& x, const Element& y) const;
        Element& axpyin(Element& r, const Element& a, const Element& x) const;

        // -- axmy:   r <- a * x - y
        // -- axmyin: r <- a * x - r
        Element& axmy  (Element& r, const Element& a, const Element& x, const Element& y) const;
        Element& axmyin(Element& r, const Element& a, const Element& x) const;

        // -- maxpy:   r <- y - a * x
        // -- maxpyin: r <- r - a * x
        Element& maxpy  (Element& r, const Element& a, const Element& x, const Element& y) const;
        Element& maxpyin(Element& r, const Element& a, const Element& x) const;

        // --------
        // -- type_string
        static const std::string type_string () {
            return "Modular<" + TypeString<Storage_t>::get()
                    + (std::is_same<Storage_t, Compute_t>::value ?
+                        "" : ", " + TypeString<Compute_t>::get() ) +  ">";
        }

        // ----- Random generators
        typedef ModularRandIter<Self_t> RandIter;
        typedef GeneralRingNonZeroRandIter<Self_t> NonZeroRandIter;
        template< class Random > Element& random(Random& g, Element& r) const
        { return init(r, g()); }
        template< class Random > Element& nonzerorandom(Random& g, Element& a) const
        {
            while (Self_t::isZero(init(a, g())));
            return a;
        }

    };
}

#include "givaro/modular-inttype.inl"

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
