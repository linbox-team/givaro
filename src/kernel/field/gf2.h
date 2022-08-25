// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Original authors (LinBox): B. Hovinen, JG Dumas, C. Pernet
// Imported and adapted by: A. Breust
// ==========================================================================

#ifndef __Givaro_field_gf2_H
#define __Givaro_field_gf2_H

#include "givaro-config.h"

#include <iostream>
#include <climits>
#include <cmath>
#include <vector>

#include "givaro/givranditer.h"
#include "givaro/udl.h"
#include "givaro/givinteger.h"

namespace Givaro
{
	/**
	 * \brief Integers modulo 2
	 *
	 * This is a tuned implementation of the field of integers modulo
	 * 2. In particular, when one constructs a VectorDomain object over
	 * this field, highly optimized bit operations will be used to make
	 * vector arithmetic very fast.
	 \ingroup field
	 */

	class GF2
	{
	public:
		const bool zero = false;
		const bool one  = true;
		const bool mOne = true;

		/** Element type
		*/
		using Self_t = GF2;
		using Residu_t = uint8_t;
		using Element = bool;
		using BitVector = std::vector<bool>;
		using BitReference = BitVector::reference;

		/** Random
		*/
		using RandIter = GIV_randIter<Self_t, bool>;
		BitReference random(const GivRandom& g, BitReference e, const Residu_t& size=0) const
		    { return e = static_cast<bool>(g() & 1u); }

		BitReference nonzerorandom(const GivRandom& g, BitReference e, const Residu_t& size=0) const
		    { return e = true; }

        Element& random(const GivRandom& g, Element& e, const Residu_t& size=0) const
		    { return e = static_cast<bool>(g() & 1u); }

		Element& nonzerorandom(const GivRandom& g, Element&  e, const Residu_t& size=0) const
		    { return e = true; }

		/** @name Object Management
		*/
		//@{

		/** Default constructor.
		*/
		GF2 () {}
		GF2 (int p, int exp = 1)
		{
			assert(p == 2);
			assert(exp == 1);
		}

		/** Copy constructor.
		 * Constructs Givaro::GF2 object by copying the field.
		 * This is required to allow field objects to be passed by value
		 * into functions.
		 * @param  F Givaro::GF2 object.
		 */
		GF2 (const GF2& F) {}

		/** Assignment operator.
		 * Required by the archetype
		 *
		 * @param F constant reference to Givaro::Modular object
		 * @return reference to Givaro::Modular object for self
		 */
		GF2& operator=(const GF2& F)
		{
			return *this;
		}

        /** Accessors
         */
		inline Element minElement() const { return false; }
		inline Element maxElement() const { return true; }

        /** Access to the modulus
         */
        inline Residu_t residu() const { return 2_ui8; }
        inline Residu_t size() const { return 2_ui8; }
        inline Residu_t characteristic() const { return 2_ui8; }
        inline Residu_t cardinality() const { return 2_ui8; }
        template<class T> inline T& characteristic(T& p) const { return p = 2u; }
        template<class T> inline T& cardinality(T& p) const { return p = 2u; }

		/** Initialization of field base element from an Integer.
		 * Behaves like C++ allocator construct.
		 * This function assumes the output field base element x has already been
		 * constructed, but that it is not already initialized.
		 * This is not a specialization of the template function because
		 * such a specialization is not allowed inside the class declaration.
		 * @return reference to field base element.
		 * @param x field base element to contain output (reference returned).
		 * @param y Integer.
		 */
		Element& init (Element& x, const int32_t &y ) const
		{
			return x = y & 1;
		}

		Element& init (Element& x, const uint32_t &y ) const
		{
			return x = y & 1;
		}

		Element& init (Element& x, const int64_t &y ) const
		{
			return x = y & 1;
		}

		Element& init (Element& x, const uint64_t &y ) const
		{
			return x = y & 1;
		}

		Element& init (Element& x, const float &y) const
		{
			return x = static_cast<unsigned char>(y) & 1;
		}

		Element& init (Element& x, const double &y) const
		{
			return x = static_cast<unsigned char>(y) & 1;
		}

		Element& init (Element& x, const Integer& y) const
		{
			return x = y & (uint32_t)1;
		}


		Element& init(Element& x) const
		{
			return x = false;
		}

		BitReference init(BitReference x, const Integer& y = 0) const
		{
			return x = y & (uint32_t)1;
		}

		/** Conversion of field base element to a template class T.
		 * This function assumes the output field base element x has already been
		 * constructed, but that it is not already initialized.
		 * @return reference to template class T.
		 * @param x template class T to contain output (reference returned).
		 * @param y constant field base element.
		 */
		Integer& convert (Integer& x, const Element& y) const
		{
			return x = y;
		}

		BitReference convert (BitReference x, const Element& y) const
		{
			return x = y;
		}

		template<class T>
		T& convert (T& x, const Element& y) const
		{
			return x = Caster<T>(y);
		}

		/** Assignment of one field base element to another.
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element& assign (Element& x, const Element& y) const
		{
			return x = y;
		}

		BitReference assign (BitReference x, const Element& y) const
		{
			return x = y;
		}

		/** Cardinality.
		 * Return Integer representing cardinality of the domain.
		 * Returns a non-negative Integer for all domains with finite
		 * cardinality, and returns -1 to signify a domain of infinite
		 * cardinality.
		 * @return Integer representing cardinality of the domain
		 */
		Integer& cardinality (Integer& c) const
		{
			return c = 2;
		}

		/** Characteristic.
		 * Return Integer representing characteristic of the domain.
		 * Returns a positive Integer to all domains with finite characteristic,
		 * and returns 0 to signify a domain of infinite characteristic.
		 * @return Integer representing characteristic of the domain.
		 */
		Integer& characteristic (Integer& c) const
		{
			return c = 2;
		}

		//@} Object Management

		/** @name Arithmetic Operations
		 * x <- y op z; x <- op y
		 * These operations require all elements, including x, to be initialized
		 * before the operation is called.  Uninitialized field base elements will
		 * give undefined results.
		 */
		//@{

		/** Equality of two elements.
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return boolean true if equal, false if not.
		 * @param  x field base element
		 * @param  y field base element
		 */
		bool areEqual (const Element& x, const Element& y) const
		{
			return x == y;
		}

		/** Zero equality.
		 * Test if field base element is equal to zero.
		 * This function assumes the field base element has already been
		 * constructed and initialized.
		 * @return boolean true if equals zero, false if not.
		 * @param  x field base element.
		 */
		bool isZero (const Element& x) const
		{
			return !x;
		}

		/** One equality.
		 * Test if field base element is equal to one.
		 * This function assumes the field base element has already been
		 * constructed and initialized.
		 * @return boolean true if equals one, false if not.
		 * @param  x field base element.
		 */
		bool isOne (const Element& x) const
		{
			return x;
		}

		/** Invertibility.
		 * Test if field base element is invertible.
		 * This function assumes the field base element has already been
		 * constructed and initialized.
		 * @return boolean true if equals one, false if not.
		 * @param  x field base element.
		 */
		bool isUnit (const Element& x) const
		{
			return x;
		}

		/** MOne equality.
		 * Test if field base element is equal to one.
		 * This function assumes the field base element has already been
		 * constructed and initialized.
		 * @return boolean true if equals one, false if not.
		 * @param  x field base element.
		 */
		bool isMOne (const Element& x) const
		{
			return x;
		}

		//@} Arithmetic Operations

		/** @name Input/Output Operations */
		//@{

		/** Print field.
		 * @return output stream to which field is written.
		 * @param  os  output stream to which field is written.
		 */
		std::ostream& write (std::ostream& os) const;

		/** Read field.
		 * @return input stream from which field is read.
		 * @param  is  input stream from which field is read.
		 */
		std::istream& read (std::istream& is);

		/** Print field base element.
		 * This function assumes the field base element has already been
		 * constructed and initialized.
		 * @return output stream to which field base element is written.
		 * @param  os  output stream to which field base element is written.
		 * @param  x   field base element.
		 */
		std::ostream& write (std::ostream& os, const Element& x) const;

		/** Read field base element.
		 * @pre This function assumes the field base element has already been
		 * constructed and initialized.
		 * @return input stream from which field base element is read.
		 * @param  is  input stream from which field base element is read.
		 * @param  x   field base element.
		 */
		std::istream& read (std::istream& is, Element& x) const;
		std::istream &read (std::istream &is, BitReference x) const;

		//@}

		/** @name Arithmetic Operations
		 * x <- y op z; x <- op y
		 * These operations require all elements, including x, to be initialized
		 * before the operation is called.  Uninitialized field base elements will
		 * give undefined results.
		 */
		//@{

		/** Addition.
		 * x = y + z
		 * This function assumes all the field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 * @param  z field base element.
		 */
		Element& add (Element& x, const Element& y, const Element& z) const;
		BitReference add (BitReference x, const Element& y, const Element& z) const;

		/** Subtraction.
		 * x = y - z
		 * This function assumes all the field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 * @param  z field base element.
		 */
		Element& sub (Element& x, const Element& y, const Element& z) const;
		BitReference sub (BitReference x, const Element& y, const Element& z) const;

		/** Multiplication.
		 * x = y * z
		 * This function assumes all the field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 * @param  z field base element.
		 */
		Element& mul (Element& x, const Element& y, const Element& z) const;
		BitReference mul (BitReference x, const Element& y, const Element& z) const;

		/** Division.
		 * x = y / z
		 * This function assumes all the field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 * @param  z field base element.
		 */
		Element& div (Element& x, const Element& y, const Element& z ) const;
		BitReference div (BitReference x, const Element& y, const Element& z ) const;

		/** Additive Inverse (Negation).
		 * x = - y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element& neg (Element& x, const Element& y) const;
		BitReference neg (BitReference x, const Element& y) const;

		/** Multiplicative Inverse.
		 * x = 1 / y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element& inv (Element& x, const Element& y) const;
		BitReference inv (BitReference x, const Element& y) const;

		/** Natural AXPY.
		 * @return reference to r.
		 * @param  r
		 * @param  a
		 * @param  x
		 * @param  y
		 */
		BitReference axpy (BitReference r, const Element& a, const Element& x, const Element& y) const;
		Element& axpy (Element& r, const Element& a, const Element& x, const Element& y) const;

		/** Natural AXMY.
		 * @return reference to r.
		 * @param  r
		 * @param  a
		 * @param  x
		 * @param  y
		 */
		BitReference axmy (BitReference r, const Element& a, const Element& x, const Element& y) const;
		Element& axmy (Element& r, const Element& a, const Element& x, const Element& y) const;

		/** Natural MAXPY.
		 * @return reference to r.
		 * @param  r
		 * @param  a
		 * @param  x
		 * @param  y
		 */
		BitReference maxpy (BitReference r, const Element& a, const Element& x, const Element& y) const;
		Element& maxpy (Element& r, const Element& a, const Element& x, const Element& y) const;

		//@} Arithmetic Operations

		/** @name Inplace Arithmetic Operations
		 * x <- x op y; x <- op x
		 */
		//@{

		/** Inplace Addition.
		 * x += y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element& addin (Element& x, const Element& y) const;
		BitReference addin (BitReference x, const Element& y) const;

		/** Inplace Subtraction.
		 * x -= y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element& subin (Element& x, const Element& y) const;
		BitReference subin (BitReference x, const Element& y) const;

		/** Inplace Multiplication.
		 * x *= y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element& mulin (Element& x, const Element& y) const;
		BitReference mulin (BitReference x, const Element& y) const;

		/** Inplace Division.
		 * x /= y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element& divin (Element& x, const Element& y) const;
		BitReference divin (BitReference x, const Element& y) const;

		/** Inplace Additive Inverse (Inplace Negation).
		 * x = - x
		 * This function assumes the field base element has already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 */
		Element& negin (Element& x) const;
		BitReference negin (BitReference x) const;

		/** Inplace Multiplicative Inverse.
		 * x = 1 / x
		 * This function assumes the field base elementhas already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 */
		Element& invin (Element& x) const;
		BitReference invin (BitReference x) const;

		/** Inplace AXPY.
		 * @return reference to r.
		 * @param  r field element (reference returned).
		 * @param  a field element.
		 * @param  x field element.
		 */
		Element& axpyin (Element& r, const Element& a, const Element& x) const;
		BitReference axpyin (BitReference r, const Element& a, const Element& x) const;

		/** Inplace AXMY.
		 * @return reference to r.
		 * @param  r field element (reference returned).
		 * @param  a field element.
		 * @param  x field element.
		 */
		Element& axmyin (Element& r, const Element& a, const Element& x) const;
		BitReference axmyin (BitReference r, const Element& a, const Element& x) const;

		/** Inplace MAXPY.
		 * @return reference to r.
		 * @param  r field element (reference returned).
		 * @param  a field element.
		 * @param  x field element.
		 */
		Element& maxpyin (Element& r, const Element& a, const Element& x) const;
		BitReference maxpyin (BitReference r, const Element& a, const Element& x) const;

		//@} Inplace Arithmetic Operations

		static inline int maxCardinality()
		{
			return 2;
		}

	}; // class GF2
}

#include "givaro/gf2.inl"

#endif


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
