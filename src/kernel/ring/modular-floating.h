// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: Clement Pernet <clement.pernet@gmail.com>
//          Pascal Giorgi  <pascal.giorgi@ens-lyon.fr>
//          A. Breust (taken from FFLAS-FFPACK)
// ==========================================================================

/*! @file field/modular-floating.h
 * @ingroup field
 * @brief  representation of <code>Z/mZ</code> over \c float.
 */

#ifndef __GIVARO_modular_floating_H
#define __GIVARO_modular_floating_H

#include <float.h>

#include "givaro/givinteger.h"
#include "givaro/givcaster.h"
#include "givaro/givranditer.h"
#include "givaro/givtypestring.h"
#include "givaro/ring-interface.h"
#include "givaro/modular-general.h"
#include "givaro/modular-implem.h"

namespace Givaro
{

    template<typename _Storage_t, typename _Compute_t>
    class Modular<_Storage_t, _Compute_t, typename std::enable_if<std::is_floating_point<_Storage_t>::value>::type>:
    public Modular_implem<_Storage_t, _Compute_t, typename make_unsigned_int<_Storage_t>::type>
    {
    public:

        using Storage_t = _Storage_t;
        using Compute_t = _Compute_t;
        using Residu_t = typename make_unsigned_int<_Storage_t>::type;

        using Element = Storage_t;
        using Self_t = Modular<_Storage_t, _Compute_t>;
        using Parent_t = Modular_implem<Storage_t, Compute_t, Residu_t>;

        // ----- Constructors
        using Modular_implem<Storage_t, Compute_t, Residu_t>::Modular_implem;
        virtual ~Modular() {}

        using Parent_t::_p;
        using Parent_t::_pc;

        inline Compute_t fcharacteristic() const { return _pc; }

        // ----- Initialisation
        Element& init (Element& x) const;

        __GIVARO_CONDITIONAL_TEMPLATE(Source,
                                      std::is_same<Source,double>::value &&
                                      std::is_same<Storage_t, float>::value)
        inline Element& init (Element&, const Source) const;

        __GIVARO_CONDITIONAL_TEMPLATE(Source,
                                      std::is_integral<Source>::value && std::is_signed<Source>::value &&
                                      sizeof(Source) >= sizeof(Storage_t))
        inline Element& init (Element&, const Source) const;

        __GIVARO_CONDITIONAL_TEMPLATE(Source,
                                      std::is_integral<Source>::value && std::is_unsigned<Source>::value &&
                                      sizeof(Source) >= sizeof(Storage_t))
        inline Element& init (Element&, const Source) const;

        inline Element& init (Element&, const Integer&) const final;

        __GIVARO_CONDITIONAL_TEMPLATE(Source,
                                      !(std::is_integral<Source>::value && sizeof(Source) >= sizeof(Storage_t)) &&
                                      !(std::is_same<Source, double>::value && std::is_same<Storage_t, float>::value) &&
                                      !std::is_same<Source, Integer&>::value)
        inline Element& init(Element& r, const Source& a) const
        { r = Caster<Element>(a); return reduce(r); }

        // ----- Reduce

        Element& reduce (Element& x, const Element& y) const;
        Element& reduce (Element& x) const;

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

        // -- type_string
        static const std::string type_string () {
            return "Modular<" + TypeString<Storage_t>::get()
                    + (sizeof(Storage_t) == sizeof(Compute_t) ?
                        "" : ", " + TypeString<Compute_t>::get() ) +  ">";
        }

        // ----- Random generators
        typedef ModularRandIter<Self_t> RandIter;
        typedef GeneralRingNonZeroRandIter<Self_t> NonZeroRandIter;
        template< class Random > Element& random(Random& g, Element& r) const
        { return init(r, g()); }
        template< class Random > Element& nonzerorandom(Random& g, Element& a) const
        { while (this->isZero(init(a, g())))
            ;
            return a; }

    };

} // Givaro

#include "modular-floating.inl"

#endif // __GIVARO_modular_floating_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
