// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: A. Breust <alexis.breust@imag.fr>
//          B. Grenet <bruno.grenet@lirmm.fr>
// ==========================================================================

#ifndef __GIVARO_modular_ruint_H
#define __GIVARO_modular_ruint_H

#include "recint/ruint.h"
#include "givaro/givinteger.h"
#include "givaro/givtypestring.h"
#include "givaro/ring-interface.h"
#include "givaro/modular-general.h"
#include "givaro/givranditer.h"
#include "givaro/modular-implem.h"

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
              using Modular_implem<Storage_t, Compute_t, Storage_t>::Modular_implem;

              using Parent_t::_p;
              using Parent_t::_pc;


              // ----- Initialisation
              Element& init (Element& x) const
              { return x = this->zero; }

              template<typename T> Element& init(Element& r, const T& a) const
              {
                  reduce(r, Caster<Element>((a < 0)? -a : a));
                  if (a < 0) negin(r);
                  return r;
              }

              Element& init(Element& r, const Integer& a) const
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
              { RecInt::rand(r); mod_n(r, _p); return r; }
              template< class Random > Element& nonzerorandom(Random& g, Element& a) const
              { while (this->isZero(random(g, a))) { } return a; }

          };

}


#include "givaro/modular-ruint.inl"

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
