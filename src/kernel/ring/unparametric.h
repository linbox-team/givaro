// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: BB <brice.boyer@lip6.fr>
//          A. Breust (taken from FFLAS-FFPACK)
// ========================================================================


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

		Element& random (Element& x) const
		{
		    x.random();
			return x;
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
