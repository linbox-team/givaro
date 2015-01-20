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

#ifndef __GIVARO_ring_zring_H
#define __GIVARO_ring_zring_H

#include <algorithm>
#include <typeinfo>

#include "givaro/unparametric-operations.h"
#include "givaro/givranditer.h"
#include "givaro/givinteger.h"

namespace Givaro
{
    /** Class ZRing
     * Ring of integers, using the templatedElement base type
     */
	template<class _Element>
	class ZRing : public UnparametricOperations<_Element> {

	public:

		/** The field's element type.
		 * Type Element must provide a default constructor,
		 * a copy constructor, a destructor, and an assignment operator.
		 */

		typedef _Element Element;
		typedef ZRing<Element> Self_t;
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

		/** Builds this field. Assumes q=0 and e = 1.
		 * Ensures consistency with field interface.
		 */
		ZRing(long int q = 0, size_t e = 1):
				one(1),zero(0),mOne(-one)
			{}  
		//@}

		template<class T>
		ZRing (const T& ) : one(1), zero(0), mOne(-one){}

		/// construct this field as copy of F.
		ZRing (const ZRing &F) : one(F.one),zero(F.zero),mOne(F.mOne){}

		Integer &cardinality (Integer &c) const {return c = 0;}

		Integer &characteristic (Integer &c) const{return c = 0;}

		int64_t &cardinality (int64_t &c) const {return c = 0;}

		int64_t &characteristic (int64_t &c) const {return c = 0;}

		int64_t cardinality () const {return 0;}

		int64_t characteristic () const {return 0;}

		ZRing<Element> operator=(const ZRing<Element> &e) {return *this ;}

		template <typename Src>
		Element& init (Element& x, const Src& s) const { return x = static_cast<const Element&>(s); }

		Element& reduce (Element& x, const Element& y) const {return init (x,y);}
		Element& reduce (Element& x) const {return init (x,x);}

		template <typename T>
		T& convert (T &x, const Element &y) const {return x = static_cast<const T&>(y);}

		    // To ensure interface consistency
		size_t minElement() const {return 0;}
		size_t maxElement() const {return 0;}

	};

} // Givaro

#endif // __FIELD_UNPARAMETRIC_H_
