// ==========================================================================
// Copyright(c)'1994-2016 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: Bastien Vialla <bastien.vialla@lirmm.fr>
// ==========================================================================

#ifndef __GIVARO_MODULAR_EXTENDED_H
#define __GIVARO_MODULAR_EXTENDED_H

#include <type_traits>

#include "givaro/givconfig.h"

#include "givaro/givranditer.h"
#include "givaro/givtypestring.h"
#include "givaro/ring-interface.h"
#include "givaro/modular-general.h"

namespace Givaro{
    /*
     *
     * Modular double/float allowing big moduli
     * !!: RandIter does not works, use your own random
     *
     */
    template<class _Element>
    class ModularExtended
    {
    public:

        typedef _Element Element;
        typedef Element* Element_ptr ;
        typedef const Element ConstElement;
        typedef const Element* ConstElement_ptr;
        // ----- Exported Types and constantes
        typedef ModularExtended<_Element> Self_t;
        using Compute_t = _Element;
        typedef uint64_t Residu_t;
        enum { size_rep = sizeof(Residu_t) };

        using is_elt_integral = std::false_type;
        static constexpr bool is_elt_integral_v = false;
        using is_elt_floating_point = std::true_type;
        static constexpr bool is_elt_floating_point_v = true;

    private:
        // Verkampt Split
        inline void split(const Element x, Element &x_h, Element &x_l) const {
            Element c;
            if(std::is_same<Element, double>::value){
                c = (Element)((1 << 27)+1);
            }else if(std::is_same<Element, float>::value){
                c = (Element)((1 << 13)+1);
            }
            c = c*x;
            x_h = c-(c-x);
            x_l = x - x_h;
        }

        // Dekker mult, a * b = s + t
        /* Note: this function may not work correctly when compiled for the 387
         * floating point coprocessor which may store temporary on 80-bit
         * registry.
         */
        inline void mult_dekker(const Element a, const Element b, Element &s, Element &t) const{
            s = a*b;
            Element ah, al, bh, bl;
            split(a, ah, al);
            split(b, bh, bl);
            t = (al*bl-(((s-ah*bh)-al*bh)-ah*bl));
        }

    public:
        // ----- Constantes
        const Element zero = 0.0;
        const Element one = 1.0;
        const Element mOne = -1.0;

        // ----- Constructors
        ModularExtended() = default;

        template<class XXX> ModularExtended(const XXX& p)
        : zero(0.0), one(1.0), mOne((Element)p - 1.0), _p((Element)p), _invp((Element)1/(Element)_p), _negp(-_p), _lp((Residu_t)p)
        {
            assert(_p >= minCardinality());
            assert(_p <= maxCardinality());
        }

        // ----- Accessors
        inline Element minElement() const { return zero; }
        inline Element maxElement() const { return mOne; }

        // ----- Access to the modulus
        inline Residu_t residu() const { return _lp; }
        inline Residu_t size() const { return _lp; }
        inline Residu_t characteristic() const { return _lp; }
        template<class T> inline T& characteristic(T& p) const { return p = _lp; }
        inline Residu_t cardinality() const { return _lp; }
        template<class T> inline T& cardinality(T& p) const { return p = _lp; }
        static inline Residu_t maxCardinality() {
            if(std::is_same<Element, double>::value)
                return 1125899906842623; // 2^(52-2) - 1
            else if(std::is_same<Element, float>::value)
                return 2097151; // 2^(23-2) - 1
        }
        static inline Residu_t minCardinality() { return 2; }

        // ----- Checkers
        inline bool isZero(const Element& a) const { return a == zero; }
        inline bool isOne (const Element& a) const { return a == one; }
        inline bool isMOne(const Element& a) const { return a == mOne; }
        inline bool isUnit(const Element& a) const;
        inline bool areEqual(const Element& a, const Element& b) const { return a == b; }
        inline size_t length(const Element a) const { return size_rep; }

        // ----- Ring-wise operators
        inline bool operator==(const Self_t& F) const { return _p == F._p; }
        inline bool operator!=(const Self_t& F) const { return _p != F._p; }
        inline Self_t& operator=(const Self_t& F)
        {
            F.assign(const_cast<Element&>(one),  F.one);
            F.assign(const_cast<Element&>(zero), F.zero);
            F.assign(const_cast<Element&>(mOne), F.mOne);
            _p = F._p;
            _negp = F._negp;
            _invp = F._invp;
            _lp= F._lp;
            return *this;
        }



        // ----- Initialisation
        inline Element& init (Element& x) const
        {
            return x = zero;
        }

        template<typename T>
        Element& init (Element& r, T a) const{
            r = Caster<Element>(a);
            return reduce(r);
        }

        Element &assign (Element &x, const Element &y) const{
            return x = y;
        }

        // ----- Convert and reduce
        template<typename T> T& convert(T& r, const Element& a) const
        { return r = static_cast<T>(a); }

        Element& reduce (Element& x, const Element& y) const{
            x=y;
            return reduce(x);
        }

        Element& reduce (Element& x) const ;

        // ----- Classic arithmetic
        Element& mul(Element& r, const Element& a, const Element& b) const;

        Element& div(Element& r, const Element& a, const Element& b) const{
            return mulin(inv(r, a), b);
        }
        Element& add(Element& r, const Element& a, const Element& b) const {
            r = a + b;
            if(r >= _p)
                r += _negp;
            return r;
        }
        Element& sub(Element& r, const Element& a, const Element& b) const {
            r = a - b;
            if(r < 0)
                r += _p;
            return r;
        }
        Element& neg(Element& r, const Element& a) const {
            r = -a;
            if(r < 0)
                r += _p;
            return r;
        }
        Element& inv(Element& x, const Element& y) const{
            invext(x,y,_p);
            if (x<0) x += _p;
            return x;
        }

        Element& mulin(Element& r, const Element& a) const {
            return mul(r, r, a);
        }
        Element& divin(Element& r, const Element& y) const{
            Element iy;
            return mulin(r, inv(iy, y));
        }
        Element& addin(Element& r, const Element& a) const {
            return add(r, r, a);
        }
        Element& subin(Element& r, const Element& a) const {
            return sub(r, r, a);
        }
        Element& negin(Element& r) const {
            return neg(r, r);
        }
        Element& invin(Element& r) const {
            return inv(r, r);
        }

        // -- axpy:   r <- a * x + y
        // -- axpyin: r <- a * x + r
        Element& axpy  (Element& r, const Element& a, const Element& x, const Element& y) const {
            Element tmp;
            mul(tmp, a, x);
            return add(r, tmp, y);
        }
        Element& axpyin(Element& r, const Element& a, const Element& x) const {
            Element tmp(r);
            return axpy(r, a, x, tmp);
        }

        // -- axmy:   r <- a * x - y
        // -- axmyin: r <- a * x - r
        Element& axmy  (Element& r, const Element& a, const Element& x, const Element& y) const {
            Element tmp;
            mul(tmp, a, x);
            return sub(r, tmp, y);
        }
        Element& axmyin(Element& r, const Element& a, const Element& x) const {
            return axmy(r, a, x, r);
        }

        // -- maxpy:   r <- y - a * x
        // -- maxpyin: r <- r - a * x
        Element& maxpy  (Element& r, const Element& a, const Element& x, const Element& y) const {
            Element tmp;
            mul(tmp, a, x);
            return sub(r, y, tmp);
        }
        Element& maxpyin(Element& r, const Element& a, const Element& x) const {
            return maxpy(r, a, x, r);
        }

        // -- type_string
        static const std::string type_string () {
            return "ModularExtended<" + TypeString<Element>::get() +  ">";
        }

        // ----- Random generators
        typedef ModularRandIter<Self_t> RandIter;
        typedef GeneralRingNonZeroRandIter<Self_t> NonZeroRandIter;
        template< class Random > Element& random(const Random& g, Element& r) const { return init(r, g()); }
        template< class Random > Element& nonzerorandom(const Random& g, Element& a) const {
            while (isZero(init(a, g())));
            return a;
        }

        // --- IO methods
        std::istream& read (std::istream& s);
        std::ostream& write(std::ostream& s) const;
        std::istream& read (std::istream& s, Element& a) const;
        std::ostream& write(std::ostream& s, const Element& a) const;

    protected:
        Element _p = 0;
        Element _invp = 0;
        Element _negp = 0;
        Residu_t _lp = 0;

    };

} // Givaro

#include "givaro/modular-extended.inl"

#endif // __GIVARO_MODULAR_EXTENDED_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
