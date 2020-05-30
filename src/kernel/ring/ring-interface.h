// ==========================================================================
// Copyright(c)'1994-2014 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: A. Breust
// ==========================================================================

#ifndef __GIVARO_ring_interface_H
#define __GIVARO_ring_interface_H

namespace Givaro
{

    /* All Givaro rings should follow this interface.
     * For instance Modular<> or Montgomery<>.
     */
    template<class _Element>
    struct RingInterface
    {
        virtual ~RingInterface() = default;

        // ----- Typedefs
        typedef _Element Element;
        typedef Element* Element_ptr ;
        typedef const Element ConstElement;
        typedef const Element* ConstElement_ptr;

        // ----- Checkers
        virtual bool isZero(const Element& a) const = 0;
        virtual bool isOne (const Element& a) const = 0;
        virtual bool isMOne(const Element& a) const = 0;
        virtual bool isUnit(const Element& a) const = 0;
        virtual bool areEqual(const Element& a, const Element& b) const = 0;

        // ----- Empty constructor
        virtual Element& init(Element& r) const = 0;

        // ----- Assignment operator
        virtual Element& assign(Element& r, const Element& a) const = 0;

        // ----- Classic arithmetic
        virtual Element& mul(Element& r, const Element& a, const Element& b) const = 0;
        virtual Element& add(Element& r, const Element& a, const Element& b) const = 0;
        virtual Element& sub(Element& r, const Element& a, const Element& b) const = 0;
        virtual Element& neg(Element& r, const Element& a) const = 0;

        virtual Element& mulin(Element& r, const Element& a) const = 0;
        virtual Element& addin(Element& r, const Element& a) const = 0;
        virtual Element& subin(Element& r, const Element& a) const = 0;
        virtual Element& negin(Element& r) const = 0;

        // -- axpy:   r <- a * x + y
        // -- axpyin: r <- a * x + r
        virtual Element& axpy  (Element& r, const Element& a, const Element& x, const Element& y) const = 0;
        virtual Element& axpyin(Element& r, const Element& a, const Element& x) const = 0;

        // -- axmy:   r <- a * x - y
        // -- axmyin: r <- a * x - r
        virtual Element& axmy  (Element& r, const Element& a, const Element& x, const Element& y) const = 0;
        virtual Element& axmyin(Element& r, const Element& a, const Element& x) const = 0;

        // -- maxpy:   r <- y - a * x
        // -- maxpyin: r <- r - a * x
        virtual Element& maxpy  (Element& r, const Element& a, const Element& x, const Element& y) const = 0;
        virtual Element& maxpyin(Element& r, const Element& a, const Element& x) const = 0;

        // --- IO methods
        virtual std::ostream& write(std::ostream& s) const =0;
        virtual std::istream& read (std::istream& s, Element& a) const =0;
        virtual std::ostream& write(std::ostream& s, const Element& a) const =0;
    }; // class RingInterface

    template<class _Element>
    struct FieldInterface : public RingInterface<_Element>
    {
        //    virtual ~FieldInterface() noexcept(true)= default;

        // ----- Division arithmetic
        virtual _Element& div(_Element& r, const _Element& a, const _Element& b) const = 0;
        virtual _Element& inv(_Element& r, const _Element& a) const = 0;
        virtual _Element& divin(_Element& r, const _Element& a) const = 0;
        virtual _Element& invin(_Element& r) const = 0;
    }; // class FieldInterface


    template<class _Element>
    struct FiniteInterface
    {
        // ----- Accessors
        virtual _Element minElement() const = 0;
        virtual _Element maxElement() const = 0;

        // 	virtual Residu_t cardinality() const = 0;
        //  virtual Residu_t maxCardinality() const = 0;
        // 	virtual Residu_t minCardinality() const = 0;

    };



    template<class _Element>
    struct FiniteFieldInterface : public FieldInterface<_Element>, public FiniteInterface<_Element> {
        //	virtual ~FiniteFieldInterface() = default;

    };
    template<class _Element>
    struct FiniteRingInterface : public  RingInterface<_Element>, public  FiniteInterface<_Element> {
        //	virtual ~FiniteRingInterface() = default;
    };


} // namespace Givaro

#endif // __GIVARO_ring_interface_H


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
