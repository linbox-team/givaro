/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
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

#include "givaro/givinteger.h"
#include "givaro/givcaster.h"
#include "givaro/givranditer.h"
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
	class Modular_implem: public FiniteFieldInterface<_Storage_t>
	{
	public:

		using Element = _Storage_t;
		using Self_t = Modular_implem<_Storage_t, _Compute_t, _Residu_t>;
		using Storage_t = _Storage_t;
		using Compute_t = _Compute_t;
		using Residu_t = _Residu_t;

		// ----- Exported Types and constantes
		enum { size_rep = sizeof(Residu_t) };

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

		static inline Residu_t minCardinality() { return 2; }

		// -- maxCardinality
		// -- Rules: Storage_t | Compute_t | maxCardinality
		// --        ----------+-----------+---------------
		// --        (u)intN_t |  uintN_t  | 2^(N/2) - 1
		// --          intN_t  | uint2N_t  | 2^(N-1) - 1
		// --         uintN_t  | uint2N_t  |   2^(N-1) - 1 //NOTE: because of invext (should be 2^N-1)
		// --         float    |  float    |  2896    // Old values (TODO: check)
		// --         double   |  double   | 94906266 // Old values (TODO: check)
		// --         Integer  |  Integer  | 0
		// --         ruint<K> |  ruint<K> | ruint<K>::maxCardinality()
		// --         ruint<K> | ruint<K+1>| (ruint<K+1>::maxCardinality()-1).Low/2

		__GIVARO_CONDITIONAL_TEMPLATE(S = Storage_t, IS_INT(S) && (sizeof(S) == sizeof(Compute_t)))
		static Residu_t maxCardinality() {
			std::size_t k = sizeof(S);
			Residu_t repunit = ~0;
			return repunit >> 4*k; // 2^(N/2) - 1 with N = bitsize(Storage_t)
		}

		__GIVARO_CONDITIONAL_TEMPLATE(S = Storage_t, IS_SINT(S) && (2*sizeof(S) == sizeof(Compute_t)))
		static Residu_t maxCardinality() {
			Residu_t repunit = ~0;
			return repunit >> 1; // 2^(N-1)-1 with N = bitsize(Storage_t)
		}

		__GIVARO_CONDITIONAL_TEMPLATE(S = Storage_t, IS_UINT(S) && (2*sizeof(S) == sizeof(Compute_t)))
		static Residu_t maxCardinality() {
			Residu_t repunit = ~0;
			return repunit >> 1; // 2^(N-1)-1 with N = bitsize(Storage_t) // NOTE: should be 2^N-1
		}

		__GIVARO_CONDITIONAL_TEMPLATE(S = Storage_t, IS_SAME(S, float))
		static Residu_t maxCardinality() { return 2896; }

		__GIVARO_CONDITIONAL_TEMPLATE(S = Storage_t, IS_SAME(S, double))
		static Residu_t maxCardinality() { return 94906266; }

		__GIVARO_CONDITIONAL_TEMPLATE(S = Storage_t, IS_SAME(S, Integer))
		static Residu_t maxCardinality() { return -1; }

		__GIVARO_CONDITIONAL_TEMPLATE(S = Storage_t, is_same_ruint<S, Compute_t>::value)
		static Residu_t maxCardinality()
		{
			return S::maxModulus();
		}

		__GIVARO_CONDITIONAL_TEMPLATE(S = Storage_t, is_smaller_ruint<S, Compute_t>::value)
		static Residu_t maxCardinality()
		{
			Residu_t max;
			return RecInt::fill_with_1(max) >> 1;
		}

		__GIVARO_CONDITIONAL_TEMPLATE(S = Storage_t, !IS_INT(S) && !IS_FLOAT(S) && !IS_SAME(S, Integer) && !is_ruint<S>::value)
		static Residu_t maxCardinality() {
			return -1;
		}

		// ----- Checkers
		inline bool isZero(const Element& a) const override { return a == zero; }
		inline bool isOne (const Element& a) const override { return a == one; }
		inline bool isMOne(const Element& a) const override { return a == mOne; }
		inline bool areEqual(const Element& a, const Element& b) const override { return a == b; }
		inline bool isUnit(const Element& a) const override
		{
			Element u,d;
			invext(u,d,a,Caster<Element>(_p));
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

		Element& assign (Element& x, const Element& y) const override
		{
			return x = y;
		}

		// ----- Convert
		template<typename T> T& convert(T& r, const Element& a) const
		{ return Caster<T,Element>(r,a); }

    	// --------
    	// -- type_string

		__GIVARO_CONDITIONAL_TEMPLATE(S = Storage_t, IS_SINT(S))
		std::string type_string() const {
			std::size_t k = sizeof(S);
			return "Modular<int" + std::to_string(8*k) + "_t>";
		}

		__GIVARO_CONDITIONAL_TEMPLATE(S = Storage_t, IS_UINT(S))
		std::string type_string() const {
			std::size_t k = sizeof(S);
			return "Modular<uint" + std::to_string(8*k) + "_t>";
		}

		__GIVARO_CONDITIONAL_TEMPLATE(S = Storage_t, IS_SAME(S, float))
		std::string type_string() const {
			return "Modular<float>";
		}

		__GIVARO_CONDITIONAL_TEMPLATE(S = Storage_t, IS_SAME(S, double))
		std::string type_string() const {
			return "Modular<double>";
		}

		__GIVARO_CONDITIONAL_TEMPLATE(S = Storage_t, IS_SAME(S, Integer))
		std::string type_string() const {
			return "Modular<Integer>";
		}

		__GIVARO_CONDITIONAL_TEMPLATE(S = Storage_t, is_ruint<S>::value)
		std::string type_string() const {
			return "Modular<RecInt::ruint<" + std::to_string(RecInt_K<S>::value) + ">>";
		}

		__GIVARO_CONDITIONAL_TEMPLATE(S = Storage_t, !IS_INT(S) && !IS_FLOAT(S) && !IS_SAME(S, Integer) && !is_ruint<S>::value)
		std::string type_string() const {
			return "<Modular<IntType>>";
		}

    	// --------
    	// ----- IO


		using FiniteFieldInterface<_Storage_t>::init;
		virtual Element& init (Element&, const Integer&) const = 0;

		inline std::ostream& write (std::ostream& s, const Element& a) const override
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
			return s << Caster<Residu_t>(a);
		}

		std::ostream& write (std::ostream& s) const override
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

		inline std::istream& read (std::istream&, Element&) const override;

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
