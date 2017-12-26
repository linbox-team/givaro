// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: A. Breust <alexis.breust@imag.fr>
//          B. Grenet <bruno.grenet@lirmm.fr>
// ==========================================================================

#ifndef __GIVARO_mod_ruint_H
#define __GIVARO_mod_ruint_H

#include "recint/ruint.h"
#include "givaro/givinteger.h"
#include "givaro/ring-interface.h"
#include "givaro/mod-general.h"
#include "givaro/givranditer.h"
#include "givaro/mod-implem.h"

namespace Givaro
{

    //! @brief The standard arithmetic in modular rings using fixed size precision.

    template<typename _Storage_t, typename _Compute_t>
    class Modular<_Storage_t, _Compute_t, 
      typename std::enable_if<is_same_ruint<_Storage_t, _Compute_t>::value 
                           || is_smaller_ruint<_Storage_t, _Compute_t>::value>::type>:
      public Modular_implem<_Storage_t, _Compute_t, _Storage_t>
    {
    public:

        // ----- Exported Types and constantes
        using Storage_t = _Storage_t;
        using Residu_t = _Storage_t;
        using Compute_t = _Compute_t;

        using Element = Storage_t;
        using Self_t = Modular<Storage_t, Compute_t>;
        using Parent_t = Modular_implem<Storage_t, Compute_t, Storage_t>;


        // ----- Constructors
        using Parent_t::Modular_implem;
        
        using Parent_t::_p;
        using Parent_t::_pc;
        //using Parent_t::_bitsizep;


        // ----- Initialisation
        Element& init (Element& x) const
        { return x = this->zero; }

        template<typename T> Element& init(Element& r, const T& a) const
        {
            reduce(r, Caster<Element>((a < 0)? -a : a));
	    if (a < 0) negin(r);
            return r;
        }

        // ----- Convert and reduce
        Element& reduce (Element& x, const Element& y) const
        { x = y % _p; return x; }
        Element& reduce (Element& x) const
        { x %= _p; return x; }

        // ----- Classic arithmetic
        Element& mul(Element& r, const Element& a, const Element& b) const override;
        Element& div(Element& r, const Element& a, const Element& b) const override;
        Element& add(Element& r, const Element& a, const Element& b) const override;
        Element& sub(Element& r, const Element& a, const Element& b) const override;
        Element& neg(Element& r, const Element& a) const override;
        Element& inv(Element& r, const Element& a) const override;

        Element& mulin(Element& r, const Element& a) const override;
        Element& divin(Element& r, const Element& a) const override;
        Element& addin(Element& r, const Element& a) const override;
        Element& subin(Element& r, const Element& a) const override;
        Element& negin(Element& r) const override;
        Element& invin(Element& r) const override;

        // -- axpy:   r <- a * x + y
        // -- axpyin: r <- a * x + r
        Element& axpy  (Element& r, const Element& a, const Element& x, const Element& y) const override;
        Element& axpyin(Element& r, const Element& a, const Element& x) const override;

        // -- axmy:   r <- a * x - y
        // -- axmyin: r <- a * x - r
        Element& axmy  (Element& r, const Element& a, const Element& x, const Element& y) const override;
        Element& axmyin(Element& r, const Element& a, const Element& x) const override;

        // -- maxpy:   r <- y - a * x
        // -- maxpyin: r <- r - a * x
        Element& maxpy  (Element& r, const Element& a, const Element& x, const Element& y) const override;
        Element& maxpyin(Element& r, const Element& a, const Element& x) const override;

        // ----- Random generators
        typedef ModularRandIter<Self_t> RandIter;
        typedef GeneralRingNonZeroRandIter<Self_t> NonZeroRandIter;
        template< class Random > Element& random(Random& g, Element& r) const
        { RecInt::rand(r); mod_n(r, _p); return r; }
        template< class Random > Element& nonzerorandom(Random& g, Element& a) const
        { while (isZero(random(g, a))) { } return a; }

    };

}


#include "givaro/mod-ruint.inl"

#endif
