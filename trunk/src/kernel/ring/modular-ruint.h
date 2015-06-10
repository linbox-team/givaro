// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: A. Breust <alexis.breust@imag.fr>
// ==========================================================================

#ifndef __GIVARO_modular_ruint_H
#define __GIVARO_modular_ruint_H

#include "recint/ruint.h"
#include "givaro/givinteger.h"
#include "givaro/ring-interface.h"
#include "givaro/modular-general.h"
#include "givaro/givranditer.h"

namespace Givaro
{

    //! @brief The standard arithmetic in modular rings using fixed size precision.

    template<size_t K>
    class Modular<RecInt::ruint<K>, RecInt::ruint<K>> : public RingInterface<RecInt::ruint<K>>
    {
    public:

        // ----- Exported Types and constantes
        using Element = typename RingInterface<RecInt::ruint<K>>::Element;
        using Self_t = Modular<RecInt::ruint<K>, RecInt::ruint<K>>;
        using Residu_t = RecInt::ruint<K>;
        using Compute_t = RecInt::ruint<K>;
        enum { size_rep = sizeof(Residu_t) };

        // ----- Representation of vector of the Element
        typedef Element* Array;

        // ----- Constantes
        const Element zero;
        const Element one;
        const Element mOne;

        // ----- Constructors
        Modular()
            :  zero(0)
            , one(1)
            , mOne(-1)
            , _p(0) {}

        Modular(const Residu_t p)
            : zero(0)
            , one(1)
            , mOne(p-1)
            , _p(p)
        {
            assert(_p >= getMinModulus());
            assert(_p <= getMaxModulus());
        }

        Modular(const Self_t& F)
            : zero(F.zero), one(F.one), mOne(F.mOne), _p(F._p) {}

        // ----- Accessors
        inline Element minElement() const override { return zero; }
        inline Element maxElement() const override { return mOne; }

        // ----- Access to the modulus
        inline Residu_t residu() const { return _p; }
        inline Residu_t size() const { return _p; }
        inline Residu_t characteristic() const { return _p; }
        inline Residu_t cardinality() const { return _p; }
        template<class T> inline T& characteristic(T& p) const { return p = _p; }
        template<class T> inline T& cardinality(T& p) const { return p = _p; }

        //std::enable_if<Compute_t>
        static inline Residu_t getMaxModulus();
        static inline Residu_t getMinModulus() { return 2; }

        // ----- Checkers
        inline bool isZero(const Element& a) const override { return a == zero; }
        inline bool isOne (const Element& a) const override { return a == one; }
        inline bool isMOne(const Element& a) const override { return a == mOne; }
        inline bool areEqual(const Element& a, const Element& b) const override { return a == b; }
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
            return r;
        }

        Element& assign (Element& x, const Element& y) const
        { return x = y; }
    
        // ----- Convert and reduce
        template<typename T> T& convert(T& r, const Element& a) const
        { return r = static_cast<T>(a); }

        Element& reduce (Element& x, const Element& y) const
        { x = y % _p; return x; }
        Element& reduce (Element& x) const
        { x %= _p; return x; }

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
        template< class Random > Element& random(const Random& g, Element& r) const     { RecInt::rand(r); mod_n(r, _p); return r; }
        template< class Random > Element& nonzerorandom(const Random& g, Element& a) const
        { while (isZero(random(g, a)));
            return a; }

        // --- IO methods
        std::istream& read (std::istream& s);
        std::ostream& write(std::ostream& s) const;
        std::istream& read (std::istream& s, Element& a) const;
        std::ostream& write(std::ostream& s, const Element a) const;

    protected:
        // -- data representation of the domain:
        Residu_t _p;
    };

} // namespace Givaro

#include "givaro/modular-ruint.inl"

#endif
