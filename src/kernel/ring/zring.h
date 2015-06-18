// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: W. J. Turner <wjturner@acm.org>
//          Bradford Hovinen <hovinen@cis.udel.edu>
//          Clement Pernet <clement.pernet@gmail.com> (inserted into FFLAS-FFPACK)
//          A. Breust (taken from FFLAS-FFPACK)
// ==========================================================================


/*! @file field/zring.h
 * @ingroup field
 * @brief  representation of a field of characteristic 0.
 */

#ifndef __GIVARO_ring_zring_H
#define __GIVARO_ring_zring_H

#include <algorithm>
#include <typeinfo>
#include <math.h>

#include "givaro/unparametric-operations.h"
#include "givaro/givranditer.h"
#include "givaro/givinteger.h"

namespace Givaro
{    
  
    /** Class ZRing.
     *  Ring of integers, using the _Element base type.
     */
    template<class _Element>
    class ZRing : public UnparametricOperations<_Element>
    {
    public:

        // ----- Exported Types and constantes
        using Element = _Element;
        using Rep = _Element;
        using Self_t = ZRing<Element>;
        using Residu_t = Element;
        using Element_ptr = Element*;
        using ConstElement_ptr = const Element*;
        enum { size_rep = sizeof(Residu_t) };
        
        const Element one  = 1;
        const Element zero = 0;
        const Element mOne = -1;

        //----- Constructors
        ZRing() {}
        ZRing(const ZRing& F) {}
        // Needed in FFLAS, when ZRing is used as delayed field.
        template<class T> ZRing(const T&) {}

        //----- Access
        Residu_t residu() const { return static_cast<Residu_t>(0); }
        Residu_t size() const { return static_cast<Residu_t>(0); }
        Residu_t cardinality() const { return static_cast<Residu_t>(0); }
        Residu_t characteristic() const { return static_cast<Residu_t>(0); }
        template<typename T> T& cardinality(T& c) const { return c = static_cast<T>(0); }
        template<typename T> T& characteristic(T& c) const { return c = static_cast<T>(0); }
        
        static inline Residu_t maxCardinality() { return -1; }
        static inline Residu_t minCardinality() { return 2; }

        //----- Ring-wise operations
        inline bool operator==(const Self_t& F) const { return true; }
        inline bool operator!=(const Self_t& F) const { return false; }
        inline ZRing<Element>& operator=(const ZRing<Element>&) { return *this; }
        // Ring tests
        bool isZero(const Element& a) const { return a == zero; }
        bool isOne (const Element& a) const { return a == one; }
        bool isMOne(const Element& a) const { return a == mOne; }

        //----- Initialisation
        Element& init(Element& x) const { return x; }
        template <typename T> Element& init(Element& x, const T& s) const
        { return x = static_cast<const Element&>(s); }
        
        Element& assign(Element& x, const Element& y) const { return x = y; }

        //----- Convert
        template <typename T> T& convert(T& x, const Element& y) const
        { return x = static_cast<const T&>(y); }
        
        Element& reduce (Element& x, const Element& y) const { return x = y; }
        Element& reduce (Element& x) const { return x; }

        // To ensure interface consistency
        size_t minElement() const { return 0; }
        size_t maxElement() const { return 0; }


        // ----- Random generators
        typedef GeneralRingRandIter<Self_t> RandIter;
        typedef GeneralRingNonZeroRandIter<Self_t> NonZeroRandIter;
        template< class Random > Element& random(const Random& g, Element& r) const
        { return init(r, g()); }
        template< class Random > Element& nonzerorandom(const Random& g, Element& a) const
        { while (isZero(init(a, g())))
                ;
            return a; }

        //----- IO
        std::ostream& write(std::ostream &os) const
        {
            return os << "ZRing<" << typeid(Element).name() << ')';
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
    
    typedef ZRing<float> FloatDomain;
    typedef ZRing<double> DoubleDomain;
}

#endif
