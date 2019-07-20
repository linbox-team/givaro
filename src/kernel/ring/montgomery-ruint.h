// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: A. Breust <alexis.breust@imag.fr>
// ==========================================================================

#ifndef __GIVARO_montgomery_ruint_H
#define __GIVARO_montgomery_ruint_H

#include "recint/ruint.h"
#include "recint/rmgmodule.h"

#include "givaro/givcaster.h"
#include "givaro/givinteger.h"
#include "givaro/givranditer.h"
#include "givaro/ring-interface.h"

namespace Givaro
{
    template<class TYPE> class Montgomery;

    //! @brief The recint-based Montgomery ring.
    //! Only odd moduli allowed
    //! An integer (a mod p) is stored as
    //! (a * r mod 2^{2^K}) with (r = 2^{2^K} mod p).

    template<size_t K>
    class Montgomery<RecInt::ruint<K>>
    {
    public:

        // ----- Exported Types and constantes
        using Element = RecInt::ruint<K>;
        using LargeElement = RecInt::ruint<K+1>;
        using Self_t = Montgomery<RecInt::ruint<K>>;
        using Residu_t = RecInt::ruint<K>;
        enum { size_rep = sizeof(Residu_t) };

        // ----- Representation of vector of the Element
        typedef Element* Array;
        using Element_ptr = Element*;
        using ConstElement_ptr = const Element*;

        // ----- Constantes
        const Element zero;
        const Element one;
        const Element mOne;

        // ----- Constructors
        Montgomery()
        :  zero(0), one(1), mOne(-1)
           , _p(0), _p1(0), _r(0), _r2(0), _r3(0)
        {}

        Montgomery(const Residu_t& p)
        : zero(0)
          , _p(p)
        {
            RecInt::arazi_qi(_p1, -_p); // p1 = -inv(p) mod 2^(2^K)
            RecInt::mod_n(_r, -_p, _p); // r = 2^(2^K) mod p

            LargeElement ltmp;
            RecInt::lmul(ltmp, _r, _r);
            RecInt::mod_n(_r2, ltmp, _p);   // r2 = r^2 mod p
            RecInt::lmul(ltmp, _r2, _r);
            RecInt::mod_n(_r3, ltmp, _p);   // r2 = r^2 mod p

            RecInt::copy(const_cast<Element&>(one), _r);
            to_mg(const_cast<Element&>(mOne), _p - 1u);

            assert( (_p & 1u) != 0u);
            assert(_p >= minCardinality());
            assert(_p <= maxCardinality());
        }

        Montgomery(const Self_t& F)
        : zero(F.zero), one(F.one), mOne(F.mOne)
          , _p(F._p), _p1(F._p1), _r(F._r), _r2(F._r2), _r3(F._r3)
        {}

        // ----- Accessors
        inline Element minElement() const { return zero; }
        inline Element maxElement() const { return mOne; }

        // ----- Access to the modulus
        inline Residu_t residu() const { return _p; }
        inline Residu_t size() const { return _p; }
        inline Residu_t characteristic() const { return _p; }
        inline Residu_t cardinality() const { return _p; }
        template<class T> inline T& characteristic(T& p) const { return p = _p; }
        template<class T> inline T& cardinality(T& p) const { return p = _p; }

        static inline Residu_t maxCardinality() { return -1; }
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
            return *this;
        }

        // ----- Initialisation
        Element& init (Element& x) const
        { return x = 0; }
        template<typename T> Element& init(Element& r, const T& a) const
        {
            reduce(r, Caster<Element>((a < 0)? -a : a));
            if (a < 0) negin(r);
            return to_mg(r);
        }
        Element& init(Element& r, const Integer& a) const
        {
            reduce(r, Caster<Element>((a < 0)? -a : a));
            if (a < 0) negin(r);
            return to_mg(r);
        }

        Element& assign (Element& x, const Element& y) const
        { return x = y; }

        // ----- Convert and reduce
        template<typename T> T& convert(T& r, const Element& a) const
        { Element tmp; return Caster<T>(r, mg_reduc(tmp, a)); }

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

        // ----- Random generators
        typedef ModularRandIter<Self_t> RandIter;
        typedef GeneralRingNonZeroRandIter<Self_t> NonZeroRandIter;
        template< class Random > Element& random(Random& g, Element& r) const
        { RecInt::rand(r); mod_n(r, _p); return r; }
        template< class Random > Element& nonzerorandom(Random& g, Element& a) const
        { while (isZero(random(g, a))) {} return a; }

        // --- IO methods
        std::istream& read (std::istream& s);
        std::ostream& write(std::ostream& s) const;
        std::istream& read (std::istream& s, Element& a) const;
        std::ostream& write(std::ostream& s, const Element& a) const;

    protected:

        // Internal montgomery-reduction. a <- b * r^{-1}
        Element& mg_reduc(Element& a, const Element& b) const;
        Element& mg_reduc(Element& a, const LargeElement& b) const;

        // a <- b * r
        Element& to_mg(Element& a, const Element& b) const;
        Element& to_mg(Element& a) const;

    protected:

        // p is the module (must be odd and > 1)
        RecInt::ruint<K> _p;
        // p1 = -inv(p) mod 2^(2^K)
        RecInt::ruint<K> _p1;
        // r = 2^(2^K) mod p
        RecInt::ruint<K> _r;
        // r2 = r^2 mod p - used to initialize elements
        RecInt::ruint<K> _r2;
        // r3 = r^3 mod p - used to compute inverse
        RecInt::ruint<K> _r3;
    };

} // namespace Givaro

#include "givaro/montgomery-ruint.inl"

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
