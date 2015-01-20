// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givzpz32std.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givzpz32std.inl,v 1.20 2011-02-02 17:16:43 bboyer Exp $
// ==========================================================================

#ifndef __GIVARO_zpz32std_INL
#define __GIVARO_zpz32std_INL

#include <cmath>

#include "modular-defines.h"

namespace Givaro {

	// ------------------------
	// ----- Classic arithmetic
	
	inline Modular<int32_t>::Element& Modular<int32_t>::mul
		(Element& r, const Element& a, const Element& b) const
	{
		return  __GIVARO_MODULAR_INTEGER_MUL(r,(int32_t)_p,(int32_t)a,(int32_t)b);
	}

	inline Modular<int32_t>::Element& Modular<int32_t>::sub
		(Element& r, const Element& a, const Element& b) const
	{
		return __GIVARO_MODULAR_INTEGER_SUB(r,(int32_t)_p,(int32_t)a,(int32_t)b);
	}

	inline Modular<int32_t>::Element& Modular<int32_t>::add
		(Element& r, const Element& a, const Element& b) const
	{
		__GIVARO_MODULAR_INTEGER_ADD(r,(int32_t)_p,(int32_t)a,(int32_t)b);
		return r;
	}

	inline Modular<int32_t>::Element& Modular<int32_t>::neg
		(Element& r, const Element& a) const
	{
		return __GIVARO_MODULAR_INTEGER_NEG(r,(int32_t)_p,(int32_t)a);
	}

	inline Modular<int32_t>::Element& Modular<int32_t>::inv
		(Element& r, const Element& a) const
	{
		int32_t u;
		Modular<int32_t>::invext(u, a, (int32_t)_p);
		return r = (u<0)?(Modular<int32_t>::Element)u + (int32_t)_p:(Modular<int32_t>::Element)u;
	}

	inline Modular<int32_t>::Element& Modular<int32_t>::div
		(Element& r, const Element& a, const Element& b) const
	{
		return mulin( inv(r,b), a );
	}
	inline Modular<int32_t>::Element& Modular<int32_t>::mulin
		(Element& r, const Element& a) const
	{
		return __GIVARO_MODULAR_INTEGER_MULIN(r,(int32_t)_p, (int32_t)a);
	}

	inline Modular<int32_t>::Element& Modular<int32_t>::divin
		(Element& r, const Element& a) const
	{
		Modular<int32_t>::Element ia;
		inv(ia, a);
		return mulin(r, ia);
	}

	inline Modular<int32_t>::Element& Modular<int32_t>::addin
		(Element& r, const Element& a) const
	{
		int32_t tmp = (int32_t)r;
		__GIVARO_MODULAR_INTEGER_ADDIN(tmp,(int32_t)_p, (int32_t)a);
		return r = Element(tmp);
	}

	inline Modular<int32_t>::Element& Modular<int32_t>::subin
		(Element& r, const Element& a) const
	{
		int32_t tmp = (int32_t)r;
		__GIVARO_MODULAR_INTEGER_SUBIN(tmp,(int32_t)_p, (int32_t)a);
		return r = (Modular<int32_t>::Element)tmp;
	}

	inline Modular<int32_t>::Element& Modular<int32_t>::negin
		(Element& r) const
	{
		return __GIVARO_MODULAR_INTEGER_NEGIN(r,(int32_t)_p);
	}

	inline Modular<int32_t>::Element& Modular<int32_t>::invin
		(Element& r) const
	{
		int32_t u;
		Modular<int32_t>::invext(u, r, (int32_t)_p);
		return r = (u<0)?(Modular<int32_t>::Element)u + (int32_t)_p:(Modular<int32_t>::Element)u;
	}
	
	inline Modular<int32_t>::Element& Modular<int32_t>::axpy
		(Element& r, const Element& a, const Element& b, const Element& c) const
	{
		return __GIVARO_MODULAR_INTEGER_MULADD(r, (int32_t)_p, (int32_t)a, (int32_t)b, (int32_t)c);
	}

	inline Modular<int32_t>::Element&  Modular<int32_t>::axpyin
		(Element& r, const Element& a, const Element& b) const
	{
		return __GIVARO_MODULAR_INTEGER_MULADDIN(r, (int32_t)_p, (int32_t)a, (int32_t)b);
	}
	
	inline Modular<int32_t>::Element& Modular<int32_t>::maxpy
		(Element& r, const Element& a, const Element& b, const Element& c) const
	{
		int32_t tmp;
		__GIVARO_MODULAR_INTEGER_MUL(tmp, (int32_t)_p, (int32_t)a, (int32_t)b);
		__GIVARO_MODULAR_INTEGER_SUB(r, (int32_t)_p, (int32_t)c, tmp);
		return r;
	}

	inline Modular<int32_t>::Element&  Modular<int32_t>::axmy
		(Element& r, const Element& a, const Element& b, const Element& c) const
	{
		int32_t tmp;
		__GIVARO_MODULAR_INTEGER_MULSUB(tmp, (int32_t)_p, (int32_t)a, (int32_t)b, (int32_t)c);
		return r = (Modular<int32_t>::Element)tmp;
	}

	inline Modular<int32_t>::Element&  Modular<int32_t>::maxpyin
		(Element& r, const Element& a, const Element& b) const
	{
		__GIVARO_MODULAR_INTEGER_SUBMULIN(r, (int32_t)_p, (int32_t)a, (int32_t)b );
		return r;
	}

	inline Modular<int32_t>::Element&  Modular<int32_t>::axmyin
		(Element& r, const Element& a, const Element& b) const
	{
		return __GIVARO_MODULAR_INTEGER_MULSUB(r, (int32_t)_p, (int32_t)a, (int32_t)b, r );
	}
	
	// ----------------------------------
	// ----- Classic arithmetic on arrays

	inline void Modular<int32_t>::mul
		(const size_t sz, Array r, constArray a, constArray b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int32_t tmp;
			__GIVARO_MODULAR_INTEGER_MUL(tmp, (int32_t)_p,(int32_t)a[i], (int32_t)b[i]);
			r[i] = (Modular<int32_t>::Element)tmp;
		}
	}

	inline void Modular<int32_t>::mul
		(const size_t sz, Array r, constArray a, Element b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int32_t tmp;
			__GIVARO_MODULAR_INTEGER_MUL(tmp, (int32_t)_p, (int32_t)a[i], (int32_t)b);
			r[i] = (Modular<int32_t>::Element)tmp;
		}
	}

	inline void Modular<int32_t>::div
		(const size_t sz, Array r, constArray a, constArray b) const
	{
		for ( size_t i=sz ; --i ; ) {
			div( r[i], a[i], b[i]);
		}
	}

	inline void Modular<int32_t>::div
		(const size_t sz, Array r, constArray a, Element b) const
	{
		Modular<int32_t>::Element ib;
		inv(ib, b);
		mul(sz, r, a, ib);
	}

	inline void Modular<int32_t>::add
		(const size_t sz, Array r, constArray a, constArray b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int32_t tmp;
			__GIVARO_MODULAR_INTEGER_ADD(tmp, (int32_t)_p, (int32_t)a[i], (int32_t)b[i]);
			r[i] = (Modular<int32_t>::Element)tmp;
		}
	}

	inline void Modular<int32_t>::add
		(const size_t sz, Array r, constArray a, Element b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int32_t tmp;
			__GIVARO_MODULAR_INTEGER_ADD(tmp,(int32_t)_p, (int32_t)a[i], (int32_t)b);
			r[i] = (Modular<int32_t>::Element)tmp;
		}
	}

	inline void Modular<int32_t>::sub
		(const size_t sz, Array r, constArray a, constArray b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int32_t tmp;
			__GIVARO_MODULAR_INTEGER_SUB(tmp, (int32_t)_p, (int32_t)a[i], (int32_t)b[i]);
			r[i] = (Modular<int32_t>::Element)tmp;
		}
	}

	inline void Modular<int32_t>::sub
		(const size_t sz, Array r, constArray a, Element b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int32_t tmp;
			__GIVARO_MODULAR_INTEGER_SUB(tmp, (int32_t)_p, (int32_t)a[i], (int32_t)b);
			r[i] = (Modular<int32_t>::Element)tmp;
		}
	}

	inline void Modular<int32_t>::neg
		(const size_t sz, Array r, constArray a) const
	{
		for ( size_t i=sz ; --i ; ) {
			int32_t tmp;
			__GIVARO_MODULAR_INTEGER_NEG(tmp, (int32_t)_p, (int32_t)a[i]);
			r[i] = (Modular<int32_t>::Element)tmp;
		}
	}

	inline void Modular<int32_t>::axpy
		(const size_t sz, Array r, constArray a, constArray x, constArray y) const
	{
		for ( size_t i=sz ; --i ; ) {
			int32_t tmp;
			__GIVARO_MODULAR_INTEGER_MULADD(tmp, (int32_t)_p, (int32_t)a[i], (int32_t)x[i], (int32_t)y[i]);
			r[i] = (Modular<int32_t>::Element)tmp;
		}
	}

	inline void Modular<int32_t>::axpyin
		(const size_t sz, Array r, constArray a, constArray x) const
	{
		for ( size_t i=sz ; --i ; ) {
			int32_t tmp = (int32_t)r[i];
			__GIVARO_MODULAR_INTEGER_MULADDIN(tmp, (int32_t)_p, (int32_t)a[i], (int32_t)x[i]);
			r[i] = (Modular<int32_t>::Element)tmp;
		}
	}

	inline void Modular<int32_t>::axmy
		(const size_t sz, Array r, constArray a, constArray x, constArray y) const
	{
		for ( size_t i=sz; i--; ) {
			int32_t tmp;
			__GIVARO_MODULAR_INTEGER_MULSUB(tmp, (int32_t)_p, (int32_t)a[i], (int32_t)x[i], (int32_t)y[i]);
			r[i] = (Modular<int32_t>::Element)tmp;
		}
	}

	inline void Modular<int32_t>::maxpyin
		(const size_t sz, Array r, constArray a, constArray x) const
	{
		for ( size_t i=sz ; --i ; ) {
			int32_t tmp = (int32_t)r[i];
			__GIVARO_MODULAR_INTEGER_SUBMULIN(tmp, (int32_t)_p, (int32_t)a[i], (int32_t)x[i]);
			r[i] = (Modular<int32_t>::Element)tmp;
		}
	}
	
	// --------------------
	// ----- Initialisation

	inline Modular<int32_t>::Element& Modular<int32_t>::init ( Element& r, const double a ) const
	{
		int sign; double ua;
		if (a < 0.0) { sign =-1; ua = -a;}
		else { ua = a; sign =1; }
		if ( ua > Signed_Trait<uint32_t>::max()){
			//     ua -= (double)floor(ua * _invdp)*_dp;
			ua = fmod(ua,_dp);
			r = (Element) ua;
		} else
			r = (Element)((ua >=_p) ? (uint32_t) ua % (uint32_t)_p : (uint32_t) ua);
		if (r && (sign ==-1)) r = (int32_t)_p - r;
		return r;
	}

	inline  Modular<int32_t>::Element&  Modular<int32_t>::init ( Element& r, const float a ) const
	{
		return init(r, (double)a);
	}



	inline  Modular<int32_t>::Element&  Modular<int32_t>::init ( Element& r, const unsigned long a ) const
	{
	       	return r = (Element)( a >= (unsigned long)_p ? a % (unsigned long)_p : a);
	}

	inline  Modular<int32_t>::Element&  Modular<int32_t>::init ( Element& r, const long a ) const
	{
		int sign;
		unsigned long ua;
		if (a <0) {
			sign =-1;
			ua = (unsigned long)-a;
		}
		else {
			ua = (unsigned long) a;
		       	sign =1;
		}
		r = Element((ua >=_p) ? ua % (uint32_t)_p : ua);
		if (r && (sign ==-1))
			r = (int32_t)_p - r;
		return r;
	}

	inline Modular<int32_t>::Element&  Modular<int32_t>::init ( Element& r, const Integer& Residu ) const
	{
		long tr;
		if (Residu <0) {
			// -a = b [p]
			// a = p-b [p]
			if ( Residu <= (Integer)(-_p) )
			       	tr = long( (-Residu) % _p) ;
			else
				tr = long(-Residu);
			if (tr)
				return r = Element(_p - (unsigned long)tr);
			else
				return r = zero;
		}
		else {
			if (Residu >= (Integer)_p ) tr =   long(Residu % _p) ;
			else tr = long(Residu);
			return r = Element(tr);
		}
	}

	inline Modular<int32_t>::Element& Modular<int32_t>::init( Element& a, const int i) const
	{
		return init(a,(long)i);
	}

	inline Modular<int32_t>::Element& Modular<int32_t>::init( Element& a, const unsigned int i) const
	{
		return init(a,(unsigned long)i);
	}


	inline void Modular<int32_t>::assign
	( const size_t sz, Array r, constArray a ) const
	{
		for ( size_t i=sz ; --i ; )
			r[i] = a[i];
	}

	inline  Modular<int32_t>::Element&  Modular<int32_t>::assign ( Element& r, const Element a ) const
	{
		return r=a;
	}

	inline void Modular<int32_t>::init
	( const size_t sz, Array r, constArray a ) const
	{
		for ( size_t i=sz ; --i ; )
			r[i] = a[i];
	}

	inline Modular<int32_t>::Element& Modular<int32_t>::init ( Element& r ) const
	{
		return r = zero;
	}

	inline Modular<int32_t>::Element&  Modular<int32_t>::dotprod
	( Element& r, const int bound, const size_t sz, constArray a, constArray b ) const
	{
		unsigned int stride = 1;
		if ((unsigned long)bound < GIVARO_MAXUINT32)
			//    stride = GIVARO_MAXULONG/((unsigned long)bound * (unsigned long)bound);
			//    hopefully stride is not unsigned long !?!
			stride = (unsigned int) (GIVARO_MAXULONG/((unsigned long)bound) / ((unsigned long)bound));
		unsigned long dot = (unsigned long)zero;
		if ((sz <10) && (sz <stride)) {
			for(  size_t i= sz; i--; )
				dot += (unsigned long)a[i] * (unsigned long)b[i];
			if (dot > _p)  return r = (Element)(dot % (unsigned long)_p);
			else  return r = (Element)dot;
		}
		size_t i_begin=0;
		stride &= (unsigned int)~0x1;
		if (stride ==0) {
			for(  size_t i= sz; --i; ) {
				dot += (unsigned long) a[i] * (unsigned long)b[i];
				if (dot>_p) dot %= _p;
			}
			return r = (Element)dot;
		}
		do {
			size_t min_sz = ((sz-i_begin) < stride ? (sz-i_begin) : stride);
			if ((min_sz & 0x1) !=0) {
				min_sz--; i_begin++;
				dot += (unsigned int) a++[min_sz] * (unsigned int)b++[min_sz];
			}
			if (min_sz > 1)
				for(  size_t i= min_sz; i>0; --i, --i, ++a, ++a, ++b, ++b )
				{
					dot += (unsigned int)a[0] * (unsigned int)b[0];
					dot += (unsigned int)a[1] * (unsigned int)b[1];
				}
			if (dot>_p) dot %= _p;
			i_begin += min_sz;
		} while (i_begin <sz);
		return r = (Element)dot;
	}

	inline Modular<int32_t>::Element&  Modular<int32_t>::dotprod
	( Element& r, const size_t sz, constArray a, constArray b ) const
	{
		return Modular<int32_t>::dotprod(r, (int32_t)_p, sz, a, b);
	}


	//  a -> r: int32_t to double
	inline void
	Modular<int32_t>::i2d ( const size_t sz, double* r, constArray a ) const
	{
		for (size_t i=0; i<sz; ++i) r[i] = a[i];
	}

	//  a -> r: double to int32_t
	inline void
	Modular<int32_t>::d2i ( const size_t sz, Array r, const double* a ) const
	{
		union d_2_l {
			double d;
			int32_t r[2];
		};
		//  static const double offset = 4503599627370496.0; // 2^52
		double offset = 4503599627370496.0; // 2^52
		for (size_t i=0; i<sz; ++i)
		{
			d_2_l tmp;
			// - normalization: put fractional part at the end of the representation
			tmp.d = a[i] + offset;
			r[i] = tmp.r[1];
			if (r[i] <(int32_t)_p) r[i] %= _p;
		}
		//    r[i] = (tmp.r[1] <_p ? tmp.r[1] : tmp.r[1]-_p);
		//    r[i] = (r[i] <_p ? r[i] : r[i]%_p);
		//    r[i] = (tmp.r[1] <_p ? tmp.r[1] : tmp.r[1]%_p);
	}



	// -- Input: (z, <_p>)
	inline std::istream& Modular<int32_t>::read (std::istream& s)
	{
		char ch;
		s >> std::ws >> ch;
		if (ch != '(')
			//    GivError::throw_error( GivBadFormat("Modular<int32_t>::read: syntax error: no '('"));
			std::cerr << "GivBadFormat(Modular<int32_t>::read: syntax error: no '('))" << std::endl;

		s >> std::ws >> ch;
		if (ch != 'z')
			//    GivError::throw_error( GivBadFormat("Modular<int32_t>::read: bad domain object"));
			std::cerr << "GivBadFormat(Modular<int32_t>::read: bad domain object))" << std::endl;

		s >> std::ws >> ch;
		if (ch != ',')
			//    GivError::throw_error( GivBadFormat("Modular<int32_t>::read: syntax error: no ','"));
			std::cerr << "GivBadFormat(Modular<int32_t>::read: syntax error: no ',')) " << std::endl;

		s >> std::ws >> _p;

		s >> std::ws >> ch;
		if (ch != ')')
			//    GivError::throw_error( GivBadFormat("Modular<int32_t>::read: syntax error: no ')'"));
			std::cerr << "GivBadFormat(Modular<int32_t>::read: syntax error: no ')')) " << std::endl;

		return s;
	}

	inline std::ostream& Modular<int32_t>::write (std::ostream& s ) const
	{
		return s << "Modular<int32_t> modulo " << residu();
	}

	inline std::istream& Modular<int32_t>::read (std::istream& s, Element& a) const
	{
		Integer tmp;
		s >> tmp;
		init(a, tmp);
		return s;
	}

	inline std::ostream& Modular<int32_t>::write (std::ostream& s, const Element a) const
	{
		return s << a;
	}

} // namespace Givaro

#endif // __GIVARO_zpz32std_INL
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
