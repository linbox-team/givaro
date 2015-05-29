// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: Clement Pernet <clement.pernet@gmail.com>
//          A. Breust (taken from FFLAS-FFPACK)
// ==========================================================================

/*! @file field/modular-balanced-int32.h
 * @ingroup field
 * @brief Balanced representation of <code>Z/mZ</code> over \c int32_t .
 * @warning NOT DEFINED for EVEN modulus
 */

#ifndef __GIVARO_modular_balanced_int32_H
#define __GIVARO_modular_balanced_int32_H

#include <cmath> // fmod
#include "givaro/givranditer.h"

#ifndef LINBOX_MAX_INT
#define LINBOX_MAX_INT INT32_MAX
#endif

#define NORMALISE(x) \
{ \
	if (x < mhalf_mod) return x += modulus; \
	else if (x > half_mod) return x -= modulus; \
}

#define NORMALISE_HI(x) \
{ \
			if (x > half_mod) x -= modulus; \
}

namespace Givaro {

	template<class TAG> class ModularBalanced;


	/// \ingroup field
	template <>
	class ModularBalanced<int32_t> {

	protected:

		int32_t modulus;
		int32_t half_mod;
		int32_t mhalf_mod;
		double modulusinv;

	public:

		typedef int32_t Element;
		typedef int32_t* Element_ptr;
		typedef const int32_t* ConstElement_ptr;

		const Element one  ;
		const Element zero ;
		const Element mOne ;

	public:

		static const bool balanced = true ;

		typedef ModularBalanced<Element> Self_t;
		typedef GeneralRingRandIter<Self_t> RandIter;
		typedef GeneralRingNonZeroRandIter<Self_t> NonZeroRandIter;

		//default modular field,taking 65521 as default modulus
		ModularBalanced () :
			modulus(65521)
			,one(1),zero(0),mOne(-1)
		{
			modulusinv = 1/(double)65521;
			half_mod = (65521 >> 1);
			mhalf_mod = half_mod-65520;
		}

		ModularBalanced (Element value, int32_t exp = 1) :
			modulus(value)
			,one(1),zero(0),mOne(-1)
		{
			half_mod = (modulus >> 1);
			mhalf_mod = half_mod-modulus+1;
			modulusinv = 1 / ((double) value);
		    assert(modulus >= getMinModulus());
		    assert(modulus <= getMaxModulus());
		}

		ModularBalanced (const ModularBalanced<Element>& mf) :
			modulus(mf.modulus),
			half_mod(mf.half_mod),
			mhalf_mod(mf.mhalf_mod),
			modulusinv(mf.modulusinv)
			,one(mf.one),zero(mf.zero),mOne(mf.mOne)
		{
		}

		ModularBalanced<Element> & assign(const ModularBalanced<Element> &F)
		{
			modulus = F.modulus;
			half_mod  = F.half_mod;
			mhalf_mod = F.mhalf_mod;
			// lmodulus   = F.lmodulus;
			modulusinv = F.modulusinv;
			F.assign(const_cast<Element&>(one),F.one);
			F.assign(const_cast<Element&>(zero),F.zero);
			F.assign(const_cast<Element&>(mOne),F.mOne);
			return *this;
		}

		const ModularBalanced<Element> &operator=(const ModularBalanced<Element> &F)
		{
			modulus = F.modulus;
			half_mod  = F.half_mod;
			mhalf_mod = F.mhalf_mod;
			// lmodulus   = F.lmodulus;
			modulusinv = F.modulusinv;
			F.assign(const_cast<Element&>(one),F.one);
			F.assign(const_cast<Element&>(zero),F.zero);
			F.assign(const_cast<Element&>(mOne),F.mOne);
			return *this;
		}

		template<class T> inline T& characteristic(T& p) const { return p = modulus; }
		template<class T> inline T& cardinality(T& p) const { return p = modulus; }

		inline uint32_t cardinality () const
		{
			return (uint32_t) modulus;
		}

		inline uint32_t characteristic () const
		{
			return (uint32_t) modulus;
		}

		template<class T> inline T& convert(T& x, const Element& y) const { return x = T(y); }

		inline std::ostream &write (std::ostream &os) const
		{
			return os << "ModularBalanced<int32_t> mod " << modulus;
		}

		inline std::istream &read (std::istream &is)
		{
			is >> modulus;
			half_mod = modulus/2;
			mhalf_mod = half_mod-modulus+1;
			modulusinv = 1 /((double) modulus );
			return is;
		}

		inline std::ostream &write (std::ostream &os, const Element &x) const
		{
			return os << x;
		}

		inline std::istream &read (std::istream &is, Element &x) const
		{
			Element tmp;
			is >> tmp;
			init(x,tmp);
			return is;
		}

		inline Element &init (Element & x, const double &y) const
		{
			x = (Element) fmod(y,(double)modulus);
			NORMALISE(x);
			return x;
		}

		inline Element &init (Element & x, const float &y) const
		{
			return init(x , (double) y);
		}

		template<class Element1>
		inline Element& init (Element & x, const Element1 &y) const
		{
			return reduce (x, Element(y));
		}

		inline Element& reduce (Element & x, const Element &y) const
		{
			x = (y % modulus);
			NORMALISE(x);
			return x;
		}

		inline Element& reduce (Element & x) const
		{
			x %= modulus;
			NORMALISE(x);
			return x;
		}

		inline Element& init(Element&x) const
		{
			return x = 0;
		}

		inline Element& init(Element& x, int64_t y ) const
		{
			x = Element(y % (int64_t)modulus);
			NORMALISE(x);
			return x;
		}

		inline Element& init(Element& x, uint32_t y ) const
		{
			x = Element(y % (uint32_t)modulus);
			NORMALISE_HI(x);
			return x;
		}

		inline Element& init (Element& x, uint64_t y) const
		{
			x = Element(y % (uint64_t)modulus);
			NORMALISE_HI(x);
			return x;
		}

		inline Element& assign(Element& x, const Element& y) const
		{
			return x = y;
		}

		inline bool areEqual (const Element &x, const Element &y) const
		{
			return x == y;
		}

		inline  bool isZero (const Element &x) const
		{
			return x == 0;
		}

		inline bool isOne (const Element &x) const
		{
			return x == 1;
		}

		inline bool isMOne (const Element &x) const
		{
			return x == mOne ;
		}

		inline Element &add (Element &x, const Element &y, const Element &z) const
		{
			x = y + z;
			NORMALISE(x);
			return x;
		}

		inline Element &sub (Element &x, const Element &y, const Element &z) const
		{
			x = y - z;
			NORMALISE(x);
			return x;
		}

		inline Element &mul (Element &x, const Element &y, const Element &z) const
		{
			Element q;

			q  = (Element) ((((double) y) * ((double) z)) * modulusinv);  // q could be off by (+/-) 1
			x = (Element) (y*z - q*modulus);

			NORMALISE(x);
			return x;
		}

		inline Element &div (Element &x, const Element &y, const Element &z) const
		{
			Element temp;
			return mul (x, y, inv(temp,z));
		}

		inline Element &neg (Element &x, const Element &y) const
		{
			return x = -y;
		}

		inline Element &inv (Element &x, const Element &y) const
		{
			Element d, t;
			XGCD(d, x, t, y, modulus);
			NORMALISE(x);
			return x;

		}

		inline Element &axpy (Element &r,
				      const Element &a,
				      const Element &x,
				      const Element &y) const
		{
			Element q;

			q  = (Element) (((((double) a) * ((double) x)) + (double)y) * modulusinv);  // q could be off by (+/-) 1
			r = (Element) (a * x + y - q*modulus);


			NORMALISE(r);

			return r;

		}

		inline Element &axmy (Element &r,
				      const Element &a,
				      const Element &x,
				      const Element &y) const
		{
			Element q;

			q  = (Element) (((((double) a) * ((double) x)) - (double)y) * modulusinv);  // q could be off by (+/-) 1
			r = (Element) (a * x - y - q*modulus);


			NORMALISE(r);
			return r;

		}

		inline Element &maxpy (Element &r,
				      const Element &a,
				      const Element &x,
				      const Element &y) const
		{
                    return negin(axmy(r,a,x,y));
		}

		inline Element &addin (Element &x, const Element &y) const
		{
			x += y;
			NORMALISE(x);
			return x;
		}

		inline Element &subin (Element &x, const Element &y) const
		{
			x -= y;
			NORMALISE(x);
			return x;
		}

		inline Element &mulin (Element &x, const Element &y) const
		{
			return mul(x,x,y);
		}

		inline Element &divin (Element &x, const Element &y) const
		{
			return div(x,x,y);
		}

		inline Element &negin (Element &x) const
		{
			return x = -x;
		}

		inline Element &invin (Element &x) const
		{
			return inv (x, x);
		}

		inline Element &axpyin (Element &r, const Element &a, const Element &x) const
		{
			Element q;

			q  = (Element) (((((double) a) * ((double) x)) + (double) r) * modulusinv);  // q could be off by (+/-) 1
			r = (Element) (a * x + r - q*modulus);


			NORMALISE(r);

                        return r;
		}

		inline Element &maxpyin (Element &r, const Element &a, const Element &x) const
		{
			Element q;
                        maxpy(q,a,x,r);
                        return assign(r,q);
		}

		static inline Element getMaxModulus()
		{
                    return 92681;
		}

		static  Element getMinModulus()	{return 3.0;}

		Element minElement() const
		{
			return mhalf_mod ;
		}

		Element maxElement() const
		{
			return half_mod ;
		}

	private:

		inline static void XGCD(Element& d, Element& s, Element& t, Element a, Element b)
		{
			Element  u, v, u0, v0, u1, v1, u2, v2, q, r;
			Element aneg = 0, bneg = 0;

			if (a < 0) {
				a = -a;
				aneg = 1;
			}

			if (b < 0) {
				b = -b;
				bneg = 1;
			}

			u1 = 1; v1 = 0;
			u2 = 0; v2 = 1;
			u = a; v = b;

			while (v != 0) {
				q = u / v;
				r = u % v;
				u = v;
				v = r;
				u0 = u2;
				v0 = v2;
				u2 =  u1 - q*u2;
				v2 = v1- q*v2;
				u1 = u0;
				v1 = v0;
			}

			if (aneg)
				u1 = -u1;

			if (bneg)
				v1 = -v1;

			d = u;
			s = u1;
			t = v1;
		}

	};

} // Givaro

#undef LINBOX_MAX_INT
#undef NORMALISE
#undef NORMALISE_HI

#endif // __GIVARO_modular_balanced_int32_H
