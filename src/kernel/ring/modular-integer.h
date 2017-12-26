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
		//using Parent_t::_bitsizep;
	


        //static inline Residu_t maxCardinality() { return -1; }

        // ----- Initialisation
        Element& init (Element& x) const;

        template<typename T> Element& init(Element& r, const T& a) const
        { r = Caster<Element>(a); return reduce(r); }

        // ----- Reduce
        Element& reduce (Element& x, const Element& y) const;
        Element& reduce (Element& x) const;
        
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
        { return init(r, g()); }
        template< class Random > Element& nonzerorandom(Random& g, Element& a) const
        { while (isZero(init(a, g())))
                ;
            return a; }

    };


    /* Specialisation for Modular<integer> field*/
    template <>
    class ModularRandIter<Modular<Integer, Integer> >
    {
    public:
        typedef Modular<Integer>  Ring;
        typedef Ring::Element Element;

        ModularRandIter(const Ring& R, const size_t& size = 0, const size_t& seed = 0) 
                : _ring(R)
        {
                // GivRandom will select a non-zero value, even if seed is 0
            GivRandom generator(seed);
            Givaro::Integer::seeding(generator());
        }
        Element& operator()(Element& elt)
        {
            // Create new random Elements
            Givaro::Integer::random_lessthan(elt,_ring.residu());

            return elt;
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
        const Ring& _ring;

    }; //  class ModularRandIter<Integer>
    
}// namespace Givaro

#include "givaro/modular-integer.inl"

#endif // __GIVARO_modular_integer_H
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
