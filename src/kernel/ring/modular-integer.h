// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: JG Dumas
//          A. Breust (taken from FFLAS-FFPACK)
// ==========================================================================

/*! @file givzpzInt.h
 * @ingroup zpz
 *  @brief Arithmetic on Z/pZ, with p a prime number in arbitrary precision.
 */

#ifndef __GIVARO_zpz_int_H
#define __GIVARO_zpz_int_H

#include "givaro/givbasictype.h"
#include "givaro/giverror.h"
#include "givaro/givinteger.h"
#include "givaro/givcaster.h"
#include "givaro/givranditer.h"
#include "givaro/modular-general.h"
#include "givaro/ring-interface.h"

namespace Givaro
{

    /*! @brief This class implement the standard arithmetic with Modulo Elements.
     * - The representation of an integer a in Zpz is the value a % p
     * .
     */
    template<>
    class Modular<Integer, Integer> : public virtual FiniteFieldInterface<Integer>
    {
    public:
        // ----- Exported Types and constantes
      typedef Modular<Integer> Self_t;
        typedef Integer Residu_t;                    // - type to store residue
        enum { size_rep = sizeof(Residu_t) };      // - size of the storage type

	// ----- Constantes
	const Element zero;
	const Element one;
	const Element mOne;

        // ----- Constructors
        ~Modular() noexcept {};
        
        Modular()
            : zero(static_cast<Element>(0))
            , one(static_cast<Element>(1))
            , mOne(static_cast<Element>(-1))
            , _p(static_cast<Residu_t>(0)) {}

        Modular(const Residu_t p)
            : zero(static_cast<Element>(0))
            , one(static_cast<Element>(1))
            , mOne(static_cast<Element>(p-1))
            , _p(static_cast<Residu_t>(p))
        {
            assert(_p >= minCardinality());
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
        
        static inline Residu_t maxCardinality() { return -1; }
        static inline Residu_t minCardinality() { return 2; }

        // ----- Checkers
        inline bool isZero(const Element& a) const override { return a == zero; }
        inline bool isOne (const Element& a) const override { return a == one; }
        inline bool isMOne(const Element& a) const override { return a == mOne; }
        inline bool isUnit(const Element& a) const;
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
        Element& init (Element& x) const;
        template<typename T> Element& init(Element& r, const T& a) const
        { r = Caster<Element>(a); return reduce(r); }

        Element& assign (Element& x, const Element& y) const;
    
        // ----- Convert and reduce
        template<typename T> T& convert(T& r, const Element& a) const
        { return r = static_cast<T>(a); }

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

        // --- IO methods
        std::ostream& write( std::ostream& s ) const;
        std::istream& read ( std::istream& s, Element& a ) const;
        std::ostream& write( std::ostream& s, const Element& a ) const;

    protected:
	
        Residu_t _p;
    };


    /* Specialisation for Modular<integer> field*/
    template <>
    class ModularRandIter<Modular<Integer> >
    {
    public:
        typedef Modular<Integer>  Ring;
        typedef Ring::Element Element;

        ModularRandIter(const Ring& R, const size_t& size = 0, const size_t& seed = 0) 
                : _ring(R)
        {
            unsigned long s=seed;
            if (! seed) {
                struct timeval tp;
                gettimeofday(&tp, 0) ;
                s = (unsigned long)(tp.tv_usec);
            }
            Givaro::Integer::seeding(s);
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

#endif // __GIVARO_zpz_int_H
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
