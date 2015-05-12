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

#ifndef __GIVARO_ring_unparametric_H
#define __GIVARO_ring_unparametric_H

#include <algorithm>
#include <typeinfo>
#include <cmath> // std::pow

#include "givaro/unparametric-operations.h"
#include "givaro/givranditer.h"
#include "givaro/givinteger.h"
#include "givaro/givcaster.h"

namespace Givaro
{
    /** Class UnparametricRing
     * Provides the implementation of a field/ring (of char >0 or 0) using
     * the arithmetic provided by UnparametricOperations.
     * UnparametricOperations is an empty class wrapping infix +,*,-,/,etc
     * operations into add, mul, sub, div, etc.
     * UnparametricRing contains characteristic and cardinality members and accessors.
     * This is used for instance to represent Z/pZ with NTL's implementation
     * - UnparametricRing<NTL::zz_p>
     * @warning: prefer ZZ<double>, ZZ<integer> over UnparametricRing<double>
     * and UnparametricRing<integer>
     */
	template<class _Element>
	class UnparametricRing : public UnparametricOperations<_Element> {
	protected:
		Givaro::Integer _p ; Givaro::Integer _card ;
	public:

		/** The field's element type.
		 * Type K must provide a default constructor,
		 * a copy constructor, a destructor, and an assignment operator.
		 */

		typedef typename UnparametricOperations<_Element>::Element Element;
        typedef Element FieldInt;
		typedef UnparametricRing<Element> Self_t;
		typedef GeneralRingRandIter<Self_t> RandIter;
		typedef GeneralRingNonZeroRandIter<Self_t> NonZeroRandIter;

		typedef Element* Element_ptr;
		typedef const Element* ConstElement_ptr;

		const Element one  ;
		const Element zero ;
		const Element mOne ;

		/** @name Field Object Basics.
		*/
		//@{

		/** Builds this field to have characteristic q and cardinality q<sup>e</sup>.
		 *  This constructor must be defined in a specialization.
		 */
		UnparametricRing(Givaro::Integer q = 0, long int e = 1) :
				_p(q), _card((q == 0) ? Givaro::Integer(-1) : Givaro::pow(q, e))
				,one(1),zero(0),mOne(-one)
			{}
		//@}

		template<class T>
		UnparametricRing (const T& q = 0) : _p(q), _card(q), one(1), zero(0), mOne(-one) {}


		/// construct this field as copy of F.
		UnparametricRing (const UnparametricRing &F) :
			_p(F._p), _card(F._card)
			,one(F.one),zero(F.zero),mOne(F.mOne)
		{
		}

		Integer &cardinality (Integer &c) const
		{
			return c = _card ;
		}

		Integer &characteristic (Integer &c) const
		{
			return c = (uint64_t)_p ;
		}

		uint64_t &cardinality (uint64_t &c) const
		{
			return c = _card ;
		}

		uint64_t &characteristic (uint64_t &c) const
		{
			return c = (uint64_t)_p ;
		}

		uint64_t cardinality () const
		{
			return (uint64_t)_card ;
		}

		uint64_t characteristic () const
		{
			return (uint64_t)_p ;
		}

		UnparametricRing<Element> operator=(const UnparametricRing<Element> &e)
		{
			_p = e.characteristic() ;
			_card = e.cardinality();
			return *this ;
		}

		template <typename Src>
		Element& init (Element& x, const Src& s) const
		{
			return x = (Element) s ;
			    //return x = static_cast<const Element&>(s);
		}

		Element& init (Element& x) const
		{
			return x = 0;
		}

		Element& reduce (Element& x, const Element& y) const {return init (x,y);}
		Element& reduce (Element& x) const {return init (x,x);}

		template <typename T>
		T& convert (T &x, const Element &y) const
		{
                    return Caster(x,y);
		}

		size_t minElement() const { return 0 ; }
		size_t maxElement() const { return _card-1; }
	};

} // Givaro

#endif // __FIELD_UNPARAMETRIC_H_
