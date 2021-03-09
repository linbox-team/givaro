// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
//          Clement Pernet <clement.pernet@gmail.com>
//          Brice Boyer (briceboyer) <boyer.brice@gmail.com> (modified)
//          A. Breust (taken from FFLAS-FFPACK)
// ==========================================================================

#ifndef __GIVARO_modular_balanced_float_H
#define __GIVARO_modular_balanced_float_H

#include "givaro/givinteger.h"
#include "givaro/givcaster.h"
#include "givaro/givranditer.h"
#include "givaro/givtypestring.h"
#include "givaro/ring-interface.h"
#include "givaro/modular-general.h"

namespace Givaro
{
    template<class TAG> class ModularBalanced;

    template <>
    class ModularBalanced<float>
    {
    public:

        // ----- Exported types
        using Self_t = ModularBalanced<float>;
        using Element = float;
        using Element_ptr = Element*;
        using ConstElement = const Element;
        using ConstElement_ptr = const Element*;
        using Residu_t = float;
        enum { size_rep = sizeof(Element) };

        // ----- Constantes
        const Element zero = 0.f;
        const Element one  = 1.f;
        const Element mOne = -1.f;

        // ----- Constructors
        ModularBalanced()
        : _p(0.f), _halfp(0.f), _mhalfp(0.f), _up(0u)
        {}

        ModularBalanced(float p)
        : _p(p)
          , _halfp(std::floor(_p / 2.f))
          , _mhalfp(_halfp - _p + 1.f)
          , _up(static_cast<uint32_t>(_p))
        {
            assert(_p >= minCardinality());
            assert(_p <= maxCardinality());
        }

        ModularBalanced(const Self_t& F)
        : zero(F.zero), one(F.one), mOne(F.mOne)
          , _p(F._p), _halfp(F._halfp), _mhalfp(F._mhalfp), _up(F._up)
        {}

        // ----- Accessors
        inline Element minElement() const { return _mhalfp; }
        inline Element maxElement() const { return _halfp; }

        // ----- Access to the modulus
        inline Residu_t residu() const { return _p; }
        inline Residu_t size() const { return _p; }
        inline Residu_t characteristic() const { return _p; }
        inline Residu_t cardinality() const { return _p; }
        template<class T> inline T& characteristic(T& p) const { return Caster(p,_p); }
        template<class T> inline T& cardinality(T& p) const { return Caster(p,_p); }

        static inline Residu_t maxCardinality() { return 8191.f; } // p=2^13-1 s.t. a*b+c fits in float, with abs(a,b,c) <= (p-1)/2
        static inline Residu_t minCardinality() { return 3.f; }

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
            _halfp = F._halfp;
            _mhalfp = F._mhalfp;
            _up = F._up;
            return *this;
        }

        // ----- Initialisation
        Element& init(Element& a) const;
        Element& init(Element& r, const float a) const;
        Element& init(Element& r, const double a) const;
        Element& init(Element& r, const int32_t a) const;
        Element& init(Element& r, const uint32_t a) const;
        Element& init(Element& r, const int64_t a) const;
        Element& init(Element& r, const uint64_t a) const;
        Element& init(Element& r, const Integer& a) const;
        template<typename T> Element& init(Element& r, const T& a) const
        { r = Caster<Element>(a); return reduce(r); }

        Element& assign(Element& r, const Element& a) const;

        // ----- Convert
        template<typename T> T& convert(T& r, const Element& a) const
        { return r = static_cast<T>(a); }

        Element& reduce(Element& r, const Element& a) const;
        Element& reduce(Element& r) const;

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
            return "ModularBalanced<" + TypeString<Element>::get() +  ">";
        }

        // ----- Random generators
        typedef ModularRandIter<Self_t> RandIter;
        typedef GeneralRingNonZeroRandIter<Self_t> NonZeroRandIter;
        template< class Random > Element& random(Random& g, Element& r) const
        { return init(r, g()); }
        template< class Random > Element& nonzerorandom(Random& g, Element& a) const
        { while (isZero(init(a, g())))
            ;
            return a; }

        // --- IO methods
        std::ostream& write(std::ostream& s) const;
        std::istream& read (std::istream& s, Element& a) const;
        std::ostream& write(std::ostream& s, const Element& a) const;

    protected:

        Element _p;
        Element _halfp;
        Element _mhalfp;
        uint32_t _up;
    };
}

#include "givaro/modular-balanced-float.inl"

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
