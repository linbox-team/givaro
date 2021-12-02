// ==========================================================================
// Copyright(c)'1994-2019 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: W. J. Turner <wjturner@acm.org>
//          Bradford Hovinen <hovinen@cis.udel.edu>
//          Clement Pernet <clement.pernet@gmail.com> (inserted into FFLAS-FFPACK)
//          A. Breust (taken from FFLAS-FFPACK)
//          J-G. Dumas (merged with InetegerDom)
// ==========================================================================


/*! @file field/zring.h
 * @ingroup field
 * @brief  representation of a field of characteristic 0.
 */

#ifndef __GIVARO_ring_zring_H
#define __GIVARO_ring_zring_H

#include <algorithm>
#include <math.h>

#include "givaro/unparametric-operations.h"
#include "givaro/givtypestring.h"
#include "givaro/givranditer.h"
#include "givaro/givcaster.h"

namespace Givaro
{
    template<typename Domain> struct DomainRandIter {
        typedef GeneralRingRandIter<Domain> RandIter;
    };

    /** Generic Class ZRing.
     *  Ring of integers, using the unparametric _Element base type.
     *  Provide unparametric with domain features
     *  ZRing<Element> is a sugar name for UnparametricZRing<Element>
     *  Also, there is a specialization, below, for ZRing<Integer>
     */
    template<class _Element>
    class UnparametricZRing : public UnparametricOperations<_Element>
    {
    public:

        // ----- Exported Types and constants
        using Element = _Element;
        using Rep = _Element;
        using Self_t = UnparametricZRing<Element>;
        using Parent_t = UnparametricOperations<Element>;
		// Unparametric have no residue,
        // this is used only for cardinality/characteristic which behave like integers
        // with a ZRing the element is supposed to behave like an integer
        using Residu_t = _Element;
        using Element_ptr = Element*;
        using ConstElement_ptr = const Element*;
        enum { size_rep = sizeof(Element) };

        const Element one  = 1;
        const Element zero = 0;
        const Element mOne = -1;

        //----- Constructors
        UnparametricZRing() {}
        UnparametricZRing(const UnparametricZRing& F) {}
        // Needed in FFLAS, when ZRing is used as delayed field.
        template<class T> UnparametricZRing(const T&) {}

        //----- Access
        Residu_t residu() const { return 0; }
        Residu_t size() const { return 0; }
        Residu_t cardinality() const { return 0; }
        Residu_t characteristic() const { return 0; }
        template<typename T> T& cardinality(T& c) const { return c = static_cast<T>(0); }
        template<typename T> T& characteristic(T& c) const { return c = static_cast<T>(0); }

        static inline Residu_t maxCardinality() { return -1; }
        static inline Residu_t minCardinality() { return 2; }

        //----- Ring-wise operations
        inline bool operator==(const Self_t& F) const { return true; }
        inline bool operator!=(const Self_t& F) const { return false; }
        inline UnparametricZRing<Element>& operator=(const UnparametricZRing<Element>&) { return *this; }
        // Ring tests
        bool isZero(const Element& a) const { return a == zero; }
        bool isOne (const Element& a) const { return a == one; }
        bool isMOne(const Element& a) const { return a == mOne; }
        bool isUnit(const Element& a) const { return isOne(a) || isMOne(a); }

        Element& abs(Element& x, const Element& a) const {return x= (a>0)? a: -a;}
        Element abs(const Element& a) const {return (a>0)? a: -a;}
        long compare (const Element& a, const Element& b) const {return (a>b)? 1: ((a<b)? -1 : 0);}

        //----- Initialisation
        Element& init(Element& x) const { return x; }
        template <typename T> Element& init(Element& x, const T& s) const
        { return Caster(x,s); }

        Element& assign(Element& x, const Element& y) const { return x = y; }

        //----- Convert
        template <typename T> T& convert(T& x, const Element& y) const
        { return Caster(x,y); }

        Element& reduce (Element& x, const Element& y) const { return x = y; }
        Element& reduce (Element& x) const { return x; }

        // To ensure interface consistency
        Element minElement() const { return 0; }
        Element maxElement() const { return 0; }


        // ----- Random generators
        typedef typename DomainRandIter<Self_t>::RandIter RandIter;
        typedef GeneralRingNonZeroRandIter<Self_t> NonZeroRandIter;
        template< class Random > Element& random(const Random& g, Element& r) const
        { return init(r, g()); }
        template< class Random > Element& nonzerorandom(const Random& g, Element& a) const
        { while (isZero(init(a, g())))
            ;
            return a; }

        // --------
        // -- type_string
        static const std::string type_string () {
            return "ZRing<" + TypeString<Element>::get() + '>';
        }

        //----- IO
        std::ostream& write(std::ostream &os) const
        {
            return os << type_string();
        }
        std::ostream& write(std::ostream &os, const Element& a) const
        {
            return os << a;
        }
        std::istream& read(std::istream &is, Element& a) const
        {
            return is >> a;
        }
    };

    template<typename Element>
    class ZRing : public UnparametricZRing<Element>
    {
        using Self_t = ZRing<Element>;
        using Parent_t = UnparametricZRing<Element>;
        using Parent_t::Parent_t; // inherit constructors
    };

    using FloatDomain = ZRing<float>;
    using DoubleDomain = ZRing<double>;

    }

    // ZRing<Integer>
#include "givaro/givinteger.h"
#endif
    /* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
    // vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
