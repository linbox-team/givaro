// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: JG Dumas
//          A. Breust (taken from FFLAS-FFPACK)
// ==========================================================================

#ifndef __GIVARO_modular_integer_H
#define __GIVARO_modular_integer_H

#include "givaro/givbasictype.h"
#include "givaro/giverror.h"
#include "givaro/givinteger.h"
#include "givaro/givcaster.h"
#include "givaro/givranditer.h"
#include "givaro/givtypestring.h"
#include "givaro/modular-general.h"
#include "givaro/ring-interface.h"
#include "givaro/modular-implem.h"

namespace Givaro
{

    /*! @brief This class implement the standard arithmetic with Modulo Elements.
     * - The representation of an integer a in Zpz is the value a % p
     * .
     */
    template<>
    class Modular<Integer>:
    public Modular_implem<Integer, Integer, Integer>
    {
    public:
        // ----- Exported Types and constantes

        using Storage_t = Integer;
        using Compute_t = Integer;
        using Residu_t = Integer;

        using Element = Storage_t;
        using Self_t = Modular<Storage_t, Compute_t>;
        using Parent_t = Modular_implem<Storage_t, Compute_t, Residu_t>;

        // Constructors
        using Parent_t::Modular_implem;

        using Parent_t::_p;
        using Parent_t::_pc;

        // ----- Initialisation
        Element& init (Element& x) const;

        Element& init(Element& r, const Integer& a) const final
        { r = Caster<Element>(a); return reduce(r); }
        template<typename T> Element& init(Element& r, const T& a) const
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
            return "Modular<" + TypeString<Integer>::get() +  ">";
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


    /* Specialisation for Modular<integer> ring*/
    template <>
    class ModularRandIter<Modular<Integer> >
    {
    public:
        typedef Modular<Integer>  Ring;
        typedef Integer Element;
        typedef Integer Residu_t;

        ModularRandIter(const Ring& R, const size_t& seed = 0, const Integer size = 0)
        : _size(size?size:R.cardinality()), _ring(R)
        {
            // GivRandom will select a non-zero value, even if seed is 0
            GivRandom generator(seed);
            Givaro::Integer::seeding(generator());
        }
        Element& operator()(Element& elt)
        {
            // Create new random Elements
	  Element tmp;
            Givaro::Integer::random_lessthan(tmp,_size);
	    return _ring.init(elt,tmp);
        }

        Element& random(Element& elt)
        {
            return this->operator()(elt);
        }
        Element operator()()
        {
            Element elt; return this->operator()(elt);
        }

        Element random()
        {
            return this->operator()();
        }

        const Ring& ring() const { return _ring; }

    private:
        const Residu_t _size;
        const Ring& _ring;

    }; //  class ModularRandIter<Integer>

}// namespace Givaro

#include "givaro/modular-integer.inl"

#endif // __GIVARO_modular_integer_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
