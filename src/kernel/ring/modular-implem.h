// ==========================================================================
// Copyright(c)'1994-2017 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: B. Grenet, R. Lebreton from existing files
// ==========================================================================

/*! @file ring/modular-implem.h
 * @ingroup ring
 * @brief Generic implementation of Modular
 */

#ifndef __GIVARO_modular_implem_H
#define __GIVARO_modular_implem_H

#include <type_traits>

#include "givaro/givinteger.h"
#include "givaro/givcaster.h"
#include "givaro/givranditer.h"
#include "givaro/givtypestring.h"
#include "givaro/ring-interface.h"
#include "givaro/modular-general.h"

#define __GIVARO_CONDITIONAL_TEMPLATE(T, ...) \
template<typename T, \
typename std::enable_if<(__VA_ARGS__), int>::type* = nullptr>
#define IS_INT(T) std::is_integral<T>::value
#define IS_SINT(T) std::is_integral<T>::value && std::is_signed<T>::value
#define IS_UINT(T) std::is_integral<T>::value && std::is_unsigned<T>::value
#define IS_FLOAT(T) std::is_floating_point<T>::value
#define IS_SAME(S,T) std::is_same<S, T>::value

namespace Givaro {

    /*! @brief This class implement the standard arithmetic with Modulo Elements.
     * - The representation of an integer a in Zpz is the value a % p
     * - m max is 46341
     * - p max is 46337
     * .
     */
    template<typename _Storage_t, typename _Compute_t, typename _Residu_t>
    class Modular_implem
    {
    public:

        using Element = _Storage_t;
        using Element_ptr = Element*;
        using ConstElement = const Element;
        using ConstElement_ptr = const Element*;
        using Self_t = Modular_implem<_Storage_t, _Compute_t, _Residu_t>;
        using Storage_t = _Storage_t;
        using Compute_t = _Compute_t;
        using Residu_t = _Residu_t;

        using is_elt_integral = std::is_integral<Element>;
        static constexpr bool is_elt_integral_v = is_elt_integral::value;
        using is_elt_floating_point = std::is_floating_point<Element>;
        static constexpr bool is_elt_floating_point_v = is_elt_floating_point::value;

        // ----- Exported Types and constantes
        enum { size_rep = sizeof(Element) };

        // ----- Constantes
        const Element zero;
        const Element one;
        const Element mOne;

        // ----- Constructors
        Modular_implem()
        : zero(static_cast<Element>(0))
          , one(static_cast<Element>(1))
          , mOne(static_cast<Element>(-1))
          , _p(static_cast<Residu_t>(0))
          , _pc(static_cast<Compute_t>(0))
        {}

        Modular_implem(const Residu_t p)
        : zero(static_cast<Element>(0))
          , one(static_cast<Element>(1))
          , mOne(static_cast<Element>(p-static_cast<Element>(1)))
          , _p(static_cast<Residu_t>(p))
          , _pc(static_cast<Compute_t>(p))
        {
            assert(_p >= minCardinality());
            assert(maxCardinality() < 1 || _p <= maxCardinality());
        }

        template<typename Source>
        Modular_implem(const Source& p)
        : zero(static_cast<Element>(0))
          , one(static_cast<Element>(1))
          , mOne(static_cast<Element>(p-1))
          , _p(static_cast<Residu_t>(p))
          , _pc(static_cast<Compute_t>(p))
        {
            assert(_p >= minCardinality());
            assert(maxCardinality() < 1 || _p <= maxCardinality());
        }

        Modular_implem(const Self_t& F)
        : zero(F.zero), one(F.one), mOne(F.mOne), _p(F._p), _pc(F._pc) {}

        virtual ~Modular_implem() = default;

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

        static inline Residu_t minCardinality() { return 2; }

        // -- maxCardinality
        // -- Goal: being able to store in Compute_t the result of x*y + z
        // --       when x, y and z belong to Storage_t
        // -- => Storage_t must store integers up to maxCardinality-1
        // -- => Compute_t must store integers up to p(p-1) where p = maxCardinality

        // -- Rules: Storage_t | Compute_t | maxCardinality
        // --        ----------+-----------+---------------
        // --        (u)intN_t |  uintN_t  | 2^(N/2)
        // --          intN_t  | uint2N_t  | 2^(N-1) - 1
        // --         uintN_t  | uint2N_t  | 2^N - 1 ; because 2^N can not be stored on Residu_t
        // --         float    |  float    | 4096: 2^12
        // --         double   |  double   | 94906266: floor(2^26 sqrt(2) + 1/2)
        // --         float    |  double   | 16777216: 2^24
        // --         Integer  |  Integer  | None
        // --         ruint<K> |  ruint<K> | ruint<K>::maxModulus (= 2^(2^(K-1)))
        // --         rint<K>  |  rint<K>  | ruint<K>::maxModulus (= 2^(2^(K-2)))
        // --         ruint<K> | ruint<K+1>| ruint<K>::maxCardinality/2 (= 2^(2^K-1)-1); because addition is done over ruint<K>

        __GIVARO_CONDITIONAL_TEMPLATE(S = Storage_t, IS_INT(S) && (sizeof(S) == sizeof(Compute_t)))
        static Residu_t maxCardinality() {
            std::size_t k = sizeof(S);
            // Residu_t repunit = ~0;
            // return repunit >> (k << 2);
            return (Residu_t)1 << (k << 2); // 2^(N/2) with N = bitsize(Storage_t)
        }

        __GIVARO_CONDITIONAL_TEMPLATE(S = Storage_t, IS_SINT(S) && (2*sizeof(S) == sizeof(Compute_t)))
        static Residu_t maxCardinality() {
            Residu_t repunit = ~0;
            return repunit >> 1; // 2^(N-1)-1 with N = bitsize(Storage_t)
        }

        __GIVARO_CONDITIONAL_TEMPLATE(S = Storage_t, IS_UINT(S) && (2*sizeof(S) == sizeof(Compute_t)))
        static Residu_t maxCardinality() {
            Residu_t repunit = ~0;
            return repunit; // 2^N-1 with N = bitsize(Storage_t)
        }

        __GIVARO_CONDITIONAL_TEMPLATE(S = Storage_t, IS_SAME(S, float) && IS_SAME(S, Compute_t))
        static Residu_t maxCardinality() { return 4096.f; }

        __GIVARO_CONDITIONAL_TEMPLATE(S = Storage_t, IS_SAME(S, float) && !IS_SAME(S, Compute_t))
        static Residu_t maxCardinality() { return 16777216.f; }

        __GIVARO_CONDITIONAL_TEMPLATE(S = Storage_t, IS_SAME(S, double))
        static Residu_t maxCardinality() { return 94906266; }

        __GIVARO_CONDITIONAL_TEMPLATE(S = Storage_t, IS_SAME(S, Integer))
        static Residu_t maxCardinality() { return -1; }

        __GIVARO_CONDITIONAL_TEMPLATE(S = Storage_t, is_same_ruint<S, Compute_t>::value)
        static Residu_t maxCardinality()
        {
            return S::maxModulus();
        }
        __GIVARO_CONDITIONAL_TEMPLATE(S = Storage_t, is_same_rint<S, Compute_t>::value)
        static Residu_t maxCardinality()
        {
	  //return typename S::Value::maxModulus();
	  return RecInt::ruint<RecInt_K<S>::value>::maxModulus()/2;
        }

        __GIVARO_CONDITIONAL_TEMPLATE(S = Storage_t, is_smaller_ruint<S, Compute_t>::value)
        static Residu_t maxCardinality()
        {
	  return Residu_t::maxCardinality()/2;
        }

      __GIVARO_CONDITIONAL_TEMPLATE(S = Storage_t, !IS_INT(S) && !IS_FLOAT(S) && !IS_SAME(S, Integer) && !is_ruint<S>::value && !is_rint<S>::value)
        static Residu_t maxCardinality() {
            return -1;
        }

        // ----- Checkers
        inline bool isZero(const Element& a) const { return a == zero; }
        inline bool isOne (const Element& a) const { return a == one; }
        inline bool isMOne(const Element& a) const { return a == mOne; }
        inline bool areEqual(const Element& a, const Element& b) const { return a == b; }
        inline bool isUnit(const Element& a) const
        {
            Element u,d;
            extended_euclid(u,d,a,Caster<Element>(_p));
            return isOne(d) || isMOne(d);
        }
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
            _pc = F._pc;
            return *this;
        }

        Element& assign (Element& x, const Element& y) const
        {
            return x = y;
        }

        // ----- Convert
        template<typename T> T& convert(T& r, const Element& a) const
        { return Caster<T,Element>(r,a); }

        // --------
        // -- type_string
        static const std::string type_string () {
            return "Modular_implem<" + TypeString<Storage_t>::get()
                             +  ", " + TypeString<Compute_t>::get()
                             +  ", " + TypeString<Residu_t>::get() + ">";
        }

        // --------
        // ----- IO

        // Needed for read (see below)
        // Thus it is declared "final" in current derived classes
        // 		using FiniteFieldInterface<_Storage_t>::init;
        virtual Element& init (Element&, const Integer&) const = 0;

        inline std::ostream& write (std::ostream& s, const Element& a) const
        {
            return this->write<Element>(s, a);
        }

        __GIVARO_CONDITIONAL_TEMPLATE(E = Element, !IS_FLOAT(E) && sizeof(E) == 1)
        inline std::ostream& write (std::ostream& s, const E& a) const
        {
            return s << int32_t(a);
        }

        __GIVARO_CONDITIONAL_TEMPLATE(E = Element, !IS_FLOAT(E) && sizeof(E) > 1)
        inline std::ostream& write (std::ostream& s, const E& a) const
        {
            return s << Caster<Element>(a);
        }

        __GIVARO_CONDITIONAL_TEMPLATE(E = Element, IS_FLOAT(E))
        inline std::ostream& write (std::ostream& s, const E& a) const
        {
            return s << Caster<typename make_signed_int<Storage_t>::type>(a);
        }

        std::ostream& write (std::ostream& s) const
        {
            return this->write<Element>(s);
        }

        __GIVARO_CONDITIONAL_TEMPLATE(E = Element, sizeof(E) == 1)
        std::ostream& write (std::ostream& s) const
        {
            return s << type_string() << " modulo " << (int32_t)residu();
        }

        __GIVARO_CONDITIONAL_TEMPLATE(E = Element, sizeof(E) > 1)
        std::ostream& write (std::ostream& s) const
        {
            return s << type_string() << " modulo " << residu();
        }

        // -- Input: (z, <_p>)
        inline std::istream& read (std::istream& s)
        {
            char ch;
            s >> std::ws >> ch;
            if (ch != '(')
                std::cerr << "GivBadFormat(Modular_implem::read: syntax error: no '('))" << std::endl;

            s >> std::ws >> ch;
            if (ch != 'z')
                std::cerr << "GivBadFormat(Modular_implem::read: bad domain object))" << std::endl;

            s >> std::ws >> ch;
            if (ch != ',')
                std::cerr << "GivBadFormat(Modular_implem::read: syntax error: no ',')) " << std::endl;

            s >> std::ws >> _p;

            s >> std::ws >> ch;
            if (ch != ')')
                std::cerr << "GivBadFormat(Modular_implem::read: syntax error: no ')')) " << std::endl;

            return s;
        }

        inline std::istream& read (std::istream&, Element&) const;

    protected:
        // -- data representation of the domain:
        Residu_t _p;
        Compute_t _pc;
    };


    template<typename _Storage_t, typename _Compute_t, typename _Residu_t>
    inline std::istream& Modular_implem<_Storage_t, _Compute_t, _Residu_t>::read (std::istream& s, Element& a) const
    {
        Integer tmp;
        s >> tmp;
        this->init(a, tmp);
        return s;
    }


}

#undef IS_INT
#undef IS_SINT
#undef IS_UINT
#undef IS_FLOAT
#undef IS_SAME

#endif // __GIVARO_modular_implem_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
