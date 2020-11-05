// ==========================================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: Bruno Grenet, from existing files
// ==========================================================================
//

/*! @file ring/modular-integral.h
 * @ingroup ring
 * @brief  representation of <code>Z/mZ</code> over int types.
 */

#ifndef __GIVARO_modular_integral_H
#define __GIVARO_modular_integral_H

#include "givaro/givinteger.h"
#include "givaro/givcaster.h"
#include "givaro/givranditer.h"
#include "givaro/givtypestring.h"
#include "givaro/ring-interface.h"
#include "givaro/modular-general.h"
#include "givaro/modular-implem.h"

#define IS_INT(T) std::is_integral<T>::value
#define IS_SINT(T) std::is_integral<T>::value && std::is_signed<T>::value
#define IS_UINT(T) std::is_integral<T>::value && std::is_unsigned<T>::value
#define IS_FLOAT(T) std::is_floating_point<T>::value

namespace Givaro {

    // -- Enabled for (Storage_t, Compute_t) either ((u)intN_t, (u)intN_t) or ((u)intN_t, (u)int2N_t)
    // -- Note that Compute_t is always converted to its unsigned version


    template<typename _Storage_t, typename _Compute_t>
    class Modular<_Storage_t, _Compute_t,
          typename std::enable_if<std::is_integral<_Storage_t>::value && std::is_integral<_Compute_t>::value
          && (sizeof(_Storage_t) == sizeof(_Compute_t) || 2*sizeof(_Storage_t) == sizeof(_Compute_t))>::type>:
          public Modular_implem<_Storage_t, typename std::make_unsigned<_Compute_t>::type, typename std::make_unsigned<_Storage_t>::type>
          {
          public:

              using Storage_t = _Storage_t;
              using Residu_t = typename std::make_unsigned<_Storage_t>::type;
              using Compute_t = typename std::make_unsigned<_Compute_t>::type;

              using Element = Storage_t;
              using Self_t = Modular<Storage_t, _Compute_t>;
              using Parent_t = Modular_implem<Storage_t, Compute_t, Residu_t>;

              // ----- Constructors
              using Modular_implem<Storage_t, Compute_t, Residu_t>::Modular_implem;

              using Parent_t::_p;
              using Parent_t::_pc;

              // ----- Initialisation

              Element& init (Element&) const;

              __GIVARO_CONDITIONAL_TEMPLATE(Source, IS_UINT(Source) && (sizeof(Source) > sizeof(Storage_t)))
              inline Element& init (Element&, const Source) const;

              __GIVARO_CONDITIONAL_TEMPLATE(Source, IS_SINT(Source) && (sizeof(Source) > sizeof(Storage_t)))
              inline Element& init (Element&, const Source) const;

              __GIVARO_CONDITIONAL_TEMPLATE(Source, IS_FLOAT(Source) && (sizeof(Source) >= sizeof(Storage_t)) && IS_SINT(Storage_t))
              inline Element& init (Element&, const Source) const;

              __GIVARO_CONDITIONAL_TEMPLATE(Source, IS_FLOAT(Source) && sizeof(Source) >= sizeof(Storage_t) && IS_UINT(Storage_t))
              inline Element& init (Element&, const Source) const;

              inline Element& init (Element&, const Integer&) const final;

              __GIVARO_CONDITIONAL_TEMPLATE(Source, IS_UINT(Storage_t)
                                            &&!(IS_INT(Source) && (sizeof(Source) > sizeof(Storage_t)))
                                            &&!(IS_FLOAT(Source) && (sizeof(Source) >= sizeof(Storage_t))))
              inline Element& init (Element&, const Source&) const;

              __GIVARO_CONDITIONAL_TEMPLATE(Source, IS_SINT(Storage_t)
                                            &&!(IS_INT(Source) && (sizeof(Source) > sizeof(Storage_t)))
                                            &&!(IS_FLOAT(Source) && (sizeof(Source) >= sizeof(Storage_t))))
              inline Element& init (Element&, const Source&) const;


              // ----- Reduce

              Element& reduce (Element&, const Element&) const;

              Element& reduce (Element&) const;

              // ------------------------
              // ----- Classic arithmetic

              Element& mul (Element&, const Element&, const Element&) const;
              Element& sub (Element&, const Element&, const Element&) const;
              Element& add (Element&, const Element&, const Element&) const;
              Element& neg (Element&, const Element&) const;
              Element& inv (Element&, const Element&) const;
              Element& div (Element&, const Element&, const Element&) const;

              Element& mulin (Element&, const Element&) const;
              Element& divin (Element&, const Element&) const;
              Element& addin (Element&, const Element&) const;
              Element& subin (Element&, const Element&) const;
              Element& negin (Element&) const;
              Element& invin (Element&) const;

              // Functions defined in modular-mulprecomp
              //
              // void precomp_p (Compute_t& invp) const
              // Element& mul_precomp_p (Element&, const Element&, const Element&, const Compute_t& invp) const
              //
              // void precomp_b (Compute_t& invb, const Element&) const
              // void precomp_b (Compute_t& invb, const Element&, const Compute_t& invp) const
              // Element& mul_precomp_b (Element&, const Element&, const Element&, const Compute_t& invb) const
#include"modular-mulprecomp.inl"

              // -- axpy:   r <- a * x + y
              // -- axpyin: r <- a * x + r
              // -- axmy:   r <- a * x - y
              // -- axmyin: r <- a * x - r
              // -- maxpy:   r <- y - a * x
              // -- maxpyin: r <- r - a * x

              Element& axpy (Element&, const Element&, const Element&, const Element&) const;
              Element& axpyin (Element&, const Element&, const Element&) const;
              Element& maxpy (Element&, const Element&, const Element&, const Element&) const;
              Element& axmy (Element&, const Element&, const Element&, const Element&) const;
              Element& maxpyin (Element&, const Element&, const Element&) const;
              Element& axmyin (Element&, const Element&, const Element&) const;

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
              template< class Random > Element& random(Random& g, Element& r, const Residu_t& size) const
              { return init(r, g() % size); }
              template< class Random > Element& nonzerorandom(Random& g, Element& a) const
              { while (this->isZero(init(a, g())))
                  ;
                  return a; }
              template< class Random > Element& nonzerorandom(Random& g, Element& a, const Residu_t& size) const
              { while (this->isZero(init(a, g() % size)))
                  ;
                  return a; }
          };
} // Givaro

#include "modular-integral.inl"

#undef IS_INT
#undef IS_SINT
#undef IS_UINT
#undef IS_FLOAT

#endif // __GIVARO_modular_integral_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
