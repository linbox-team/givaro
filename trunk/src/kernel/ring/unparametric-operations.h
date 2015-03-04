/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* field/unparametric.h
 * Copyright (C) 1999-2005 William J Turner,
 *               2001 Bradford Hovinen
 *               2005 Clement Pernet
 *
 * Written by W. J. Turner <wjturner@acm.org>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 * Modified By C. Pernet and inserted into Fflas_Ffpack
 *
 * ========LICENCE========
 * This file is part of the library FFLAS-FFPACK.
 *
 * FFLAS-FFPACK is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *
 */


/*! @file field/unparametric.h
 * @ingroup field
 * @brief  representation of a field of characteristic 0.
 */

#ifndef __GIVARO_ring_unparametric_operations_H
#define __GIVARO_ring_unparametric_operations_H

#include <iostream> // std::ostream
#include <string>

namespace Givaro
{

	/** \brief Unparameterized field adapter.
	 * \ingroup field
	 * \defgroup UnparametricRing UnparametricRing
	 *
	 * A field having an interface similar to that of floats is adapted to LinBox.
	 *
	 *  Used to generate efficient field classes for unparameterized fields (or hidden parameter fields).
	 *
	 *  Some fields are implemented by definition of the C++ arithmetic operators, such as z = x*y,
	 *  for z, y, z instances of a type K.   The LinBox field
	 *  Unparametric<K> is the adaptation to LinBox.
	 *
	 *  For a typical unparametric field, some of the methods must be defined in a specialization.
	 */
	template <class K>
	class UnparametricOperations {
	public:
		typedef K Element;
        
		UnparametricOperations(){}
		//@{
		virtual ~UnparametricOperations () {}

		/* Assignment operator.
		 * Assigns UnparametricRing object F to field.
		 * @param  F UnparametricRing object.
		 */
		// I believe this should be virtual -bds
		///
		//@} Field Object Basics.

		/** @name Data Object Management.
		 * first argument is set and the value is also returned.
		 */
		//@{
		Element& init (Element& x) const
		{
			return x;
		}


		///
		Element &assign (Element &x, const Element &y) const
		{
			return x = y;
		}
		//@}

		/// @name Comparison Predicates
		//@{
		///  x == y
		bool areEqual (const Element &x, const Element &y) const
		{
			return x == y;
		}

		///  x == 0
		bool isZero (const Element &x) const
		{
			return x == Element (0);
		}

		///  x == 1
		bool isOne (const Element &x) const
		{
			return x == Element (1);
		}

		inline bool isMOne (const Element &x) const
		{
			return x == Element(-1) ;
		}

		//@} Comparison Predicates


		/** @name Arithmetic Operations
		 * The first argument is set and is also the return value.
		 */
		//@{

		/// x := y + z
		Element &add (Element &x, const Element &y, const Element &z) const
		{
			return x = y + z;
		}

		/// x := y - z
		Element &sub (Element &x, const Element &y, const Element &z) const
		{
			return x = y - z;
		}

		/// x := y*z
		Element &mul (Element &x, const Element &y, const Element &z) const
		{
			return x = y * z;
		}

		/// x := y/z
		Element &div (Element &x, const Element &y, const Element &z) const
		{
			return x = y / z;
		}

		/// x := -y
		Element &neg (Element &x, const Element &y) const
		{
			return x = - y;
		}

		/// x := 1/y
		Element &inv (Element &x, const Element &y) const
		{
			return x = Element (1) / y;
		}

		/// z := a*x + y
		// more optimal implementation, if available, can be defined in a template specialization.
		Element &axpy (Element &z,
			       const Element &a,
			       const Element &x,
			       const Element &y) const
		{
			return z = a * x + y;
		}

		/// z := a*x - y
		// more optimal implementation, if available, can be defined in a template specialization.
		Element &axmy (Element &z,
			       const Element &a,
			       const Element &x,
			       const Element &y) const
		{
			return z = a * x - y;
		}

		/// z := y - a*x
		// more optimal implementation, if available, can be defined in a template specialization.
		Element &maxpy (Element &z,
			       const Element &a,
			       const Element &x,
			       const Element &y) const
		{
			return z = y - a * x;
		}

		//@} Arithmetic Operations

		/** @name Inplace Arithmetic Operations
		 * The first argument is modified and the result is the return value.
		 */
		//@{

		/// x := x + y
		Element &addin (Element &x, const Element &y) const
		{
			return x += y;
		}

		/// x := x - y
		Element &subin (Element &x, const Element &y) const
		{
			return x -= y;
		}

		/// x := x*y
		Element &mulin (Element &x, const Element &y) const
		{
			return x *= y;
		}

		/// x := x/y
		Element &divin (Element &x, const Element &y) const
		{
			return x /= y;
		}

		/// x := -x
		Element &negin (Element &x) const
		{
			return x = - x;
		}

		/// x := 1/x
		Element &invin (Element &x) const
		{
			return x = Element (1) / x;
		}

		/// y := a*x + y
		Element &axpyin (Element &y, const Element &a, const Element &x) const
		{
			return y += a * x;
		}

		/// y := a*x + y
		Element &maxpyin (Element &y,
			       const Element &a,
			       const Element &x) const
		{
			return y += a * x;
		}

		//@} Inplace Arithmetic Operations

		/** @name Input/Output Operations */
		//@{

		/** Print field.
		 * @return output stream to which field is written.
		 * @param  os  output stream to which field is written.
		 */
		std::ostream &write (std::ostream &os) const
		{
			return os << "UnparametricRing<" << sizeof(Element) <<',' << typeid(Element).name() << ')';
		}

		std::ostream &write (std::ostream &os, std::string F) const
		{
            return write(F != "" ? os << F : os);
		}

		/** Read field.
		 * @return input stream from which field is read.
		 * @param  is  input stream from which field is read.
		 */
		std::istream &read (std::istream &is) const
		{
			return is;
		}

		/** Print field element.
		 * @return output stream to which field element is written.
		 * @param  os  output stream to which field element is written.
		 * @param  x   field element.
		 */
		std::ostream &write (std::ostream &os, const Element &x) const
		{
			return os << x;
		}
		
		/** Read field element.
		 * @return input stream from which field element is read.
		 * @param  is  input stream from which field element is read.
		 * @param  x   field element.
		 */
		virtual std::istream &read (std::istream &is, Element &x) const
		{
			return is >> x;
		}

		//@}
	};

} // Givaro

#endif // __FIELD_UNPARAMETRIC_H_
