// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: JG Dumas (from P. Zimmermann's Montgomery implementation)
//          A. Breust (adapted)
// ==========================================================================

/*! @file givmontg32.h
 * @ingroup zpz
 * @brief NO DOC
 */

#ifndef __GIVARO_montg32_H
#define __GIVARO_montg32_H

#include "givaro/givbasictype.h"
#include "givaro/giverror.h"
#include "givaro/givcaster.h"
#include "givaro/givranditer.h"
#include "givaro/modular-general.h" // invext()
#include "givaro/ring-interface.h" // invext()

#include <cmath>

#define B32 65536UL
#define MASK32 65535UL
#define HALF_BITS32 16

namespace Givaro
{
    template<class TYPE> class Montgomery;

    /*! @brief This class implements the standard arithmetic with Modulo Elements.
     *   Reduction is made through Montgomery's reduction.
     *   Representation of a is by storing (aB).
     *   - We must have gcd(p,2)=1
     *   - We must have \f$(p-1)^2 + p(B-1) < B^2 \f$, i.e. \f$2<p \leq 40504\f$ for \f$B=2^16\f$.
     *   - m max is 40503
     *   - p max is 40499
     */
    template<>
    class Montgomery<int32_t> : public RingInterface<uint32_t>
    {
    public:
        // ----- Exported Types and constantes
        using Self_t = Montgomery<int32_t>;
        using Residu_t = uint32_t;
        enum { size_rep = sizeof(Residu_t) };

        // ----- Constantes
        const Element zero = 0;
        const Element one;
        const Element mOne;

        // ----- Constructors
        Montgomery() : one(1UL), mOne(0UL), _p(0UL), _dp(0.0) {}

        Montgomery( Residu_t p, int = 1) :
            one(0UL), mOne(0UL),
            _p(  (Residu_t)  p),
            _Bp( (Residu_t)  B32%p),
            _B2p((Residu_t)  (_Bp<<HALF_BITS32) % p),
            _B3p((Residu_t)  (_B2p<<HALF_BITS32) % p),
            _nim((Residu_t)  -invext(uint64_t(_p), B32)),
            _dp( (double)    p)
        {
            const_cast<Element&>(one) = _Bp;
            const_cast<Element&>(mOne) = _p - one;
        }

        Montgomery( const Montgomery<int32_t>& F)
            : one(F.one), mOne(F.mOne)
            , _p(F._p), _Bp(F._Bp), _B2p( F._B2p), _B3p( F._B3p)
            , _nim(F._nim), _dp(F._dp)
        {}

        // ----- Accessors
        inline Element minElement() const final { return zero; }
        inline Element maxElement() const final { return mOne; }

        // ----- Access to the modulus
        inline Residu_t residu() const { return _p; }
        inline Residu_t size() const { return _p; }
        inline Residu_t characteristic() const { return _p; }
        inline Residu_t cardinality() const { return _p; }
        template<class T> inline T& characteristic(T& p) const { return p = _p; }
        template<class T> inline T& cardinality(T& p) const { return p = _p; }
        static inline Residu_t getMaxModulus() { return 40503; } // 2^15.3
        static inline Residu_t getMinModulus() { return 2; }

        // ----- Checkers
        inline bool isZero(const Element& a) const final { return a == zero; }
        inline bool isOne (const Element& a) const final { return a == one; }
        inline bool isMOne(const Element& a) const final { return a == mOne; }
        inline bool areEqual(const Element& a, const Element& b) const final { return a == b; }
        inline size_t length(const Element a) const { return size_rep; }

        // ----- Ring-wise operators
        bool operator==(const Self_t& F) const { return _p == F._p; }
        bool operator!=(const Self_t& F) const { return _p != F._p; }
        Self_t& operator=(const Self_t& F)
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
        Element& init (Element& x, const double a) const;
        Element& init (Element& x, const int64_t a) const;
        Element& init (Element& x, const uint64_t a) const;
        Element& init (Element& x, const Integer& a) const;
        template<typename T> Element& init(Element& r, const T& a) const
        {
            reduce(r, Caster<Element>((a < 0)? -a : a));
	    if (a < 0) negin(r);
            return redc(r, r * _B2p);
        }

        Element& assign (Element& x, const Element& y) const
        { return x = y; }
    
        // ----- Convert and reduce
        template<typename T> T& convert(T& r, const Element& a) const
        { Element c; return r = Caster<T>(redc(c, a)); }

        Element& reduce (Element& x, const Element& y) const
        { x = y % _p; return x; }
        Element& reduce (Element& x) const
        { x %= _p; return x; }

        // ----- Classic arithmetic
        Element& mul(Element& r, const Element& a, const Element& b) const final;
        Element& div(Element& r, const Element& a, const Element& b) const final;
        Element& add(Element& r, const Element& a, const Element& b) const final;
        Element& sub(Element& r, const Element& a, const Element& b) const final;
        Element& neg(Element& r, const Element& a) const final;
        Element& inv(Element& r, const Element& a) const final;

        Element& mulin(Element& r, const Element& a) const final;
        Element& divin(Element& r, const Element& a) const final;
        Element& addin(Element& r, const Element& a) const final;
        Element& subin(Element& r, const Element& a) const final;
        Element& negin(Element& r) const final;
        Element& invin(Element& r) const final;

        // -- axpy:   r <- a * x + y
        // -- axpyin: r <- a * x + r
        Element& axpy  (Element& r, const Element& a, const Element& x, const Element& y) const final;
        Element& axpyin(Element& r, const Element& a, const Element& x) const final;

        // -- axmy:   r <- a * x - y
        // -- axmyin: r <- a * x - r
        Element& axmy  (Element& r, const Element& a, const Element& x, const Element& y) const final;
        Element& axmyin(Element& r, const Element& a, const Element& x) const final;

        // -- maxpy:   r <- y - a * x
        // -- maxpyin: r <- r - a * x
        Element& maxpy  (Element& r, const Element& a, const Element& x, const Element& y) const final;
        Element& maxpyin(Element& r, const Element& a, const Element& x) const final;

        // ----- Random generators
        typedef ModularRandIter<Self_t> RandIter;
        typedef GeneralRingNonZeroRandIter<Self_t> NonZeroRandIter;
        template< class Random > Element& random(const Random& g, Element& r) const
        { return init(r, g()); }
        template< class Random > Element& nonzerorandom(const Random& g, Element& a) const
        { while (isZero(init(a, g())));
            return a; }

        // --- IO methods
        std::ostream& write(std::ostream& s) const;
        std::istream& read (std::istream& s, Element& a) const;
        std::ostream& write(std::ostream& s, const Element& a) const;

    protected:

        Element& redc(Element&, const Element) const ;
        Element redcal(const Element) const;
        Element redcsal(const Element) const;
        Element& redcin(Element&) const;
        Element& redcs(Element&, const Element) const;
        Element& redcsin(Element&) const;

        // -- data representation of the domain:
        Residu_t _p;
        Residu_t _Bp;
        Residu_t _B2p;
        Residu_t _B3p;
        Residu_t _nim;
        double _dp;
    };
}

#include "givaro/montgomery-int32.inl"

#endif

