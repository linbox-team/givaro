// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givzpz16std.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givzpz16std.inl,v 1.18 2011-02-02 17:16:43 bboyer Exp $
// ==========================================================================
// Description:

// ---------
// -- normalized operations
// ---------
#ifndef __GIVARO_modular_int8_INL
#define __GIVARO_modular_int8_INL

#include "modular-defines.h"

namespace Givaro {

	// -------------
	// ----- Modular 

    template<>
	inline Modular<int8_t, int8_t>::Residu_t
	Modular<int8_t, int8_t>::getMaxModulus() { return 15u; } // 2^4 - 1

    template<>
	inline Modular<int8_t, uint8_t>::Residu_t
	Modular<int8_t, uint8_t>::getMaxModulus() { return 15u; }

    template<>
	inline Modular<int8_t, uint16_t>::Residu_t
	Modular<int8_t, uint16_t>::getMaxModulus() { return 127u; } // 2^7 - 1

    template<>
	inline Modular<int8_t, int16_t>::Residu_t
	Modular<int8_t, int16_t>::getMaxModulus() { return 127u; }

	// ------------------------
	// ----- Classic arithmetic
	
	template<typename COMP>
	inline typename Modular<int8_t, COMP>::Element& Modular<int8_t, COMP>::mul
		(Element& r, const Element& a, const Element& b) const
	{
		return  __GIVARO_MODULAR_INTEGER_MUL(r,_p,a,b);
	}

	template<typename COMP>
    inline typename Modular<int8_t, COMP>::Element& Modular<int8_t, COMP>::sub
		(Element& r, const Element& a, const Element& b) const
	{
		return __GIVARO_MODULAR_INTEGER_SUB(r,_p,a,b);
	}

	template<typename COMP>
    inline typename Modular<int8_t, COMP>::Element& Modular<int8_t, COMP>::add
		(Element& r, const Element& a, const Element& b) const
	{
		__GIVARO_MODULAR_INTEGER_ADD(r,_p,a,b);
		return r;
	}

	template<typename COMP>
    inline typename Modular<int8_t, COMP>::Element& Modular<int8_t, COMP>::neg
		(Element& r, const Element& a) const
	{
		return __GIVARO_MODULAR_INTEGER_NEG(r,_p,a);
	}

	template<typename COMP>
    inline typename Modular<int8_t, COMP>::Element& Modular<int8_t, COMP>::inv
		(Element& r, const Element& a) const
	{
		invext(r, a, Element(_p));
		return (r < 0)? r = Element(r + _p) : r;
	}

	template<typename COMP>
    inline typename Modular<int8_t, COMP>::Element& Modular<int8_t, COMP>::div
		(Element& r, const Element& a, const Element& b) const
	{
		return mulin(inv(r, b), a);
	}
	
	template<typename COMP>
    inline typename Modular<int8_t, COMP>::Element& Modular<int8_t, COMP>::mulin
		(Element& r, const Element& a) const
	{
		return __GIVARO_MODULAR_INTEGER_MULIN(r, _p, a);
	}

	template<typename COMP>
    inline typename Modular<int8_t, COMP>::Element& Modular<int8_t, COMP>::divin
		(Element& r, const Element& a) const
	{
		Element ia;
		return mulin(r, inv(ia, a));
	}

	template<typename COMP>
    inline typename Modular<int8_t, COMP>::Element& Modular<int8_t, COMP>::addin
		(Element& r, const Element& a) const
	{
		__GIVARO_MODULAR_INTEGER_ADDIN(r,_p,a);
		return r;
	}

	template<typename COMP>
    inline typename Modular<int8_t, COMP>::Element& Modular<int8_t, COMP>::subin
		(Element& r, const Element& a) const
	{
		__GIVARO_MODULAR_INTEGER_SUBIN(r,_p,a);
		return r;
	}

	template<typename COMP>
    inline typename Modular<int8_t, COMP>::Element& Modular<int8_t, COMP>::negin
		(Element& r) const
	{
		return __GIVARO_MODULAR_INTEGER_NEGIN(r,_p);
	}

	template<typename COMP>
    inline typename Modular<int8_t, COMP>::Element& Modular<int8_t, COMP>::invin
		(Element& r) const
	{
		return inv(r, r);
	}
	
	template<typename COMP>
    inline typename Modular<int8_t, COMP>::Element& Modular<int8_t, COMP>::axpy
		(Element& r, const Element& a, const Element& b, const Element& c) const
	{
		return __GIVARO_MODULAR_INTEGER_MULADD(r,_p,a,b,c);
	}

	template<typename COMP>
    inline typename Modular<int8_t, COMP>::Element&  Modular<int8_t, COMP>::axpyin
		(Element& r, const Element& a, const Element& b) const
	{
		return __GIVARO_MODULAR_INTEGER_MULADDIN(r,_p,a,b);
	}
	
	template<typename COMP>
    inline typename Modular<int8_t, COMP>::Element& Modular<int8_t, COMP>::maxpy
		(Element& r, const Element& a, const Element& b, const Element& c) const
	{
		Element tmp;
		__GIVARO_MODULAR_INTEGER_MUL(tmp,_p,a,b);
		__GIVARO_MODULAR_INTEGER_SUB(r,_p,c,tmp);
		return r;
	}

	template<typename COMP>
    inline typename Modular<int8_t, COMP>::Element&  Modular<int8_t, COMP>::axmy
		(Element& r, const Element& a, const Element& b, const Element& c) const
	{
		return __GIVARO_MODULAR_INTEGER_MULSUB(r,_p,a,b,c);
	}

	template<typename COMP>
    inline typename Modular<int8_t, COMP>::Element&  Modular<int8_t, COMP>::maxpyin
		(Element& r, const Element& a, const Element& b) const
	{
		__GIVARO_MODULAR_INTEGER_SUBMULIN(r,_p,a,b);
		return r;
	}

	template<typename COMP>
    inline typename Modular<int8_t, COMP>::Element&  Modular<int8_t, COMP>::axmyin
		(Element& r, const Element& a, const Element& b) const
	{
		return __GIVARO_MODULAR_INTEGER_MULSUB(r,_p,a,b,r);
	}

	// ----------------------------------
	// ----- Classic arithmetic on arrays

	template<typename COMP>
    inline void Modular<int8_t, COMP>::mul
		(const size_t sz, Array r, constArray a, constArray b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int16_t tmp;
			__GIVARO_MODULAR_INTEGER_MUL(tmp,_p,a[i],b[i]);
			r[i] = (Modular<int8_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int8_t, COMP>::mul
		(const size_t sz, Array r, constArray a, Element b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int16_t tmp;
			__GIVARO_MODULAR_INTEGER_MUL(tmp,_p,a[i],b);
			r[i] = (Modular<int8_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int8_t, COMP>::div
		(const size_t sz, Array r, constArray a, constArray b) const
	{
		for ( size_t i=sz ; --i ; ) {
			div( r[i], a[i], b[i]);
		}
	}

	template<typename COMP>
    inline void Modular<int8_t, COMP>::div
		(const size_t sz, Array r, constArray a, Element b) const
	{
		Modular<int8_t, COMP>::Element ib;
		inv(ib, b);
		mul(sz, r, a, ib);
	}

	template<typename COMP>
    inline void Modular<int8_t, COMP>::add
		(const size_t sz, Array r, constArray a, constArray b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int16_t tmp;
			__GIVARO_MODULAR_INTEGER_ADD(tmp,_p,a[i],b[i]);
			r[i] = (Modular<int8_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int8_t, COMP>::add
		(const size_t sz, Array r, constArray a, Element b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int8_t tmp;
			__GIVARO_MODULAR_INTEGER_ADD(tmp,_p,a[i],b);
			r[i] = (Modular<int8_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int8_t, COMP>::sub
		(const size_t sz, Array r, constArray a, constArray b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int8_t tmp;
			__GIVARO_MODULAR_INTEGER_SUB(tmp,_p,a[i],b[i]);
			r[i] = (Modular<int8_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int8_t, COMP>::sub
		(const size_t sz, Array r, constArray a, Element b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int8_t tmp;
			__GIVARO_MODULAR_INTEGER_SUB(tmp,_p,a[i],b);
			r[i] = (Modular<int8_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int8_t, COMP>::neg
		(const size_t sz, Array r, constArray a) const
	{
		for ( size_t i=sz ; --i ; ) {
			int8_t tmp;
			__GIVARO_MODULAR_INTEGER_NEG(tmp,_p,a[i]);
			r[i] = (Modular<int8_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int8_t, COMP>::axpy
		(const size_t sz, Array r, constArray a, constArray x, constArray y) const
	{
		for ( size_t i=sz ; --i ; ) {
			int8_t tmp;
			__GIVARO_MODULAR_INTEGER_MULADD(tmp,_p,a[i],x[i],y[i]);
			r[i] = (Modular<int8_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int8_t, COMP>::axpyin
		(const size_t sz, Array r, constArray a, constArray x) const
	{
		for ( size_t i=sz ; --i ; ) {
			int8_t tmp = (int8_t)r[i];
			__GIVARO_MODULAR_INTEGER_MULADDIN(tmp,_p,a[i],x[i]);
			r[i] = (Modular<int8_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int8_t, COMP>::axmy
		(const size_t sz, Array r, constArray a, constArray x, constArray y) const
	{
		for ( size_t i=sz; i--; ) {
			int8_t tmp;
			__GIVARO_MODULAR_INTEGER_MULSUB(tmp,_p,a[i],x[i],y[i]);
			r[i] = (Modular<int8_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int8_t, COMP>::maxpyin
		(const size_t sz, Array r, constArray a, constArray x) const
	{
		for ( size_t i=sz ; --i ; ) {
			int8_t tmp = (int8_t)r[i];
			__GIVARO_MODULAR_INTEGER_SUBMULIN(tmp,_p,a[i],x[i]);
			r[i] = (Modular<int8_t, COMP>::Element)tmp;
		}
	}
	
	// --------------------
	// ----- Initialisation
	
	template<typename COMP>
    inline typename Modular<int8_t, COMP>::Element &  Modular<int8_t, COMP>::init ( Element& r, const unsigned long a ) const
	{
		return r = (Element)( a >= (unsigned long)_p ? a % (unsigned long)_p : a);
	}

	template<typename COMP>
    inline typename Modular<int8_t, COMP>::Element &  Modular<int8_t, COMP>::init ( Element& r, const long a ) const
	{
		int sign; long ua;
		if (a <0) {
			sign =-1;
			ua = -a;
		}
		else {
			ua = a;
			sign =1;
		}
		r = Element( (ua >=_p) ? ua % (uint16_t)_p : ua );
		if (r && (sign ==-1))
			r = (Element)(_p - r);
		return r;
	}

	template<typename COMP>
    inline typename Modular<int8_t, COMP>::Element&  Modular<int8_t, COMP>::init ( Element& r, const Integer& Residu ) const
	{
		Element tr;
		if (Residu <0) {
			// -a = b [p]
			// a = p-b [p]
			if ( Residu <= (Integer)(-_p) ) tr = Element( (-Residu) % (uint16_t)_p) ;
			else tr = Element(-Residu);
			if (tr)
				return r = Element((uint16_t)_p - (uint16_t)tr);
			else
				return r = zero;
		}
		else {
			if (Residu >= (Integer)_p ) tr =   Element(Residu % _p) ;
			else tr = Element(Residu);
			return r = (Element)tr;
		}
	}




	template<typename COMP>
    inline typename Modular<int8_t, COMP>::Element & Modular<int8_t, COMP>::init( Element& a, const int i) const
	{
		return init(a,(long)i);
	}
	template<typename COMP>
    inline typename Modular<int8_t, COMP>::Element & Modular<int8_t, COMP>::init( Element& a, const double i) const
	{
		return init(a,(long)i);
	}
	template<typename COMP>
    inline typename Modular<int8_t, COMP>::Element & Modular<int8_t, COMP>::init( Element& a, const float i) const
	{
		return init(a,(double)i);
	}
	template<typename COMP>
    inline typename Modular<int8_t, COMP>::Element & Modular<int8_t, COMP>::init( Element& a, const unsigned int i) const
	{
		return init(a,(unsigned long)i);
	}


	template<typename COMP>
    inline void Modular<int8_t, COMP>::assign
	( const size_t sz, Array r, constArray a ) const
	{
		for ( size_t i=sz ; --i ; )
			r[i] = a[i];
	}

	template<typename COMP>
    inline typename Modular<int8_t, COMP>::Element &  Modular<int8_t, COMP>::assign ( Element& r, const Element a ) const
	{  return r=a;
	}

	template<typename COMP>
    inline void Modular<int8_t, COMP>::init
	( const size_t sz, Array r, constArray a ) const
	{
		for ( size_t i=sz ; --i ; )
			r[i] = a[i];
	}

	template<typename COMP>
    inline typename Modular<int8_t, COMP>::Element& Modular<int8_t, COMP>::init ( Element& r ) const
	{
		return r = zero;
	}

	template<typename COMP>
    inline typename Modular<int8_t, COMP>::Element& Modular<int8_t, COMP>::dotprod
	( Element& r, const int bound, const size_t sz, constArray a, constArray b ) const
	{
		unsigned int stride = 1;
		if ((unsigned long)bound < GIVARO_MAXUINT16)
			stride = (unsigned int) (GIVARO_MAXULONG/((unsigned long)bound * (unsigned long)bound));
		unsigned long dot = (unsigned long)zero;
		if ((sz <10) && (sz <stride)) {
			for(  size_t i= sz; i--; )
				dot += (unsigned long)a[i] * (unsigned long)b[i];
			if (dot > _p) r = (Element)(dot % (uint16_t)_p);
			else r = (Element)dot;
			return r;
		}
		unsigned int i_begin=0;
		stride &= (unsigned int)~0x1;
		if (stride ==0) {
			for(  size_t i = sz; i--; ) {
				dot += (unsigned long)a[i] * (unsigned long)b[i];
				if (dot>_p) dot %= _p;
			}
			r = (Element)dot;
			return r;
		}
		do {
			size_t min_sz = ((sz-i_begin) < stride ? (sz-i_begin) : stride);
			if ((min_sz & 0x1) !=0) {
				min_sz--;
				i_begin++;
				dot += (unsigned long)a++[min_sz] * (unsigned long)b++[min_sz];
			}
			if (min_sz > 1)
				for(  size_t i= min_sz; i>0; --i, --i, ++a, ++a, ++b, ++b )
				{
					dot += (unsigned long)a[0] * (unsigned long)b[0];
					dot += (unsigned long)a[1] * (unsigned long)b[1];
				}
			if (dot>_p) dot %= (uint16_t)_p;
			i_begin += (unsigned int) min_sz;
		} while (i_begin <sz);
		
		return r = (Element)dot;
	}

	template<typename COMP>
    inline typename Modular<int8_t, COMP>::Element& Modular<int8_t, COMP>::dotprod
	( Element& r, const size_t sz, constArray a, constArray b ) const
	{
		return Modular<int8_t, COMP>::dotprod(r, _p, sz, a, b);
	}


	//  a -> r: int16_t to double
	template<typename COMP>
	inline void Modular<int8_t, COMP>::i2d ( const size_t sz, double* r, constArray a ) const
	{
		for (size_t i=0; i<sz; ++i) r[i] = a[i];
	}

	//  a -> r: double to int16_t
	template<typename COMP>
	inline void Modular<int8_t, COMP>::d2i ( const size_t sz, Array r, const double* a ) const
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
			r[i] = (Element) tmp.r[1];
			if (r[i] <_p)
				r[i] = Element(r[i]%_p);
		}
		//    r[i] = (tmp.r[1] <_p ? tmp.r[1] : tmp.r[1]-_p);
		//    r[i] = (r[i] <_p ? r[i] : r[i]%_p);
		//    r[i] = (tmp.r[1] <_p ? tmp.r[1] : tmp.r[1]%_p);
	}



	// -- Input: (z, <_p>)
	template<typename COMP>
	inline std::istream& Modular<int8_t, COMP>::read (std::istream& s)
	{
		char ch;
		s >> std::ws >> ch;
		if (ch != '(')
			//    GivError::throw_error( GivBadFormat("Modular<int8_t, COMP>::read: syntax error: no '('"));
			std::cerr << "GivBadFormat(Modular<int8_t, COMP>::read: syntax error: no '('))" << std::endl;

		s >> std::ws >> ch;
		if (ch != 'z')
			//    GivError::throw_error( GivBadFormat("Modular<int8_t, COMP>::read: bad domain object"));
			std::cerr << "GivBadFormat(Modular<int8_t, COMP>::read: bad domain object))" << std::endl;

		s >> std::ws >> ch;
		if (ch != ',')
			//    GivError::throw_error( GivBadFormat("Modular<int8_t, COMP>::read: syntax error: no ','"));
			std::cerr << "GivBadFormat(Modular<int8_t, COMP>::read: syntax error: no ',')) " << std::endl;

		s >> std::ws >> _p;

		s >> std::ws >> ch;
		if (ch != ')')
			//    GivError::throw_error( GivBadFormat("Modular<int8_t, COMP>::read: syntax error: no ')'"));
			std::cerr << "GivBadFormat(Modular<int8_t, COMP>::read: syntax error: no ')')) " << std::endl;

		return s;
	}

	template<>
	inline std::ostream& Modular<int8_t, int8_t>::write (std::ostream& s) const
	{
		return s << "Modular<int8_t, uint8_t> modulo " << residu();
	}

	template<>
	inline std::ostream& Modular<int8_t, uint8_t>::write (std::ostream& s) const
	{
		return s << "Modular<int8_t, uint8_t> modulo " << residu();
	}

	template<>
	inline std::ostream& Modular<int8_t, int16_t>::write (std::ostream& s) const
	{
		return s << "Modular<int8_t, uint16_t> modulo " << residu();
	}

	template<>
	inline std::ostream& Modular<int8_t, uint16_t>::write (std::ostream& s) const
	{
		return s << "Modular<int8_t, uint16_t> modulo " << residu();
	}

	template<typename COMP>
	inline std::istream& Modular<int8_t, COMP>::read (std::istream& s, Element& a) const
	{
        	Integer tmp;
		s >> tmp;
		init(a, tmp);
		return s;
	}

	template<typename COMP>
	inline std::ostream& Modular<int8_t, COMP>::write (std::ostream& s, const Element a) const
	{
		// Cast needed to be understood as a number
		return s << int32_t(a);
	}

} // namespace Givaro

#endif // __GIVARO_modular_int8_INL
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
