// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: W. J. Turner <wjturner@acm.org>
//          Bradford Hovinen <hovinen@cis.udel.edu>
//          Clement Pernet <clement.pernet@gmail.com> (inserted into FFLAS-FFPACK)
//          A. Breust (taken from FFLAS-FFPACK)
// ==========================================================================


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

		uint64_t &cardinality (uint64_t &c) const {return c = 0U;}
		uint64_t &characteristic (uint64_t &c) const {return c = 0U;}
		uint64_t cardinality () const {return 0U;}
		uint64_t characteristic () const {return 0U;}

		ZRing<Element> operator=(const ZRing<Element> &e) {return *this ;}

		Element& init (Element& x) const { return x;}

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
	/* Representations of Z with floating point elements*/
	typedef ZRing<float> FloatDomain;
	typedef ZRing<double> DoubleDomain;


} // Givaro

#endif // __FIELD_UNPARAMETRIC_H_
