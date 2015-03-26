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

	// -------------
	// ----- Modular 

    template<>
	inline Modular<int32_t, int32_t>::Residu_t
	Modular<int32_t, int32_t>::getMaxModulus() { return 65535u; } // 2^16 - 1

    template<>
	inline Modular<int32_t, uint32_t>::Residu_t
	Modular<int32_t, uint32_t>::getMaxModulus() { return 65535u; }

    template<>
	inline Modular<int32_t, uint64_t>::Residu_t
	Modular<int32_t, uint64_t>::getMaxModulus() { return 2147483647u; } // 2^31 - 1

    template<>
	inline Modular<int32_t, int64_t>::Residu_t
	Modular<int32_t, int64_t>::getMaxModulus() { return 2147483647u; }

	// ------------------------
	// ----- Classic arithmetic
	
	template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::mul
		(Element& r, const Element& a, const Element& b) const
	{
		return  __GIVARO_MODULAR_INTEGER_MUL(r,_p,a,b);
	}

	template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::sub
		(Element& r, const Element& a, const Element& b) const
	{
		return __GIVARO_MODULAR_INTEGER_SUB(r,_p,a,b);
	}

	template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::add
		(Element& r, const Element& a, const Element& b) const
	{
		__GIVARO_MODULAR_INTEGER_ADD(r,_p,a,b);
		return r;
	}

	template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::neg
		(Element& r, const Element& a) const
	{
		return __GIVARO_MODULAR_INTEGER_NEG(r,_p,a);
	}

	template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::inv
		(Element& r, const Element& a) const
	{
		invext(r, a, int32_t(_p));
		return (r < 0)? r += int32_t(_p) : r;
	}

	template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::div
		(Element& r, const Element& a, const Element& b) const
	{
		return mulin(inv(r, b), a);
	}
	template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::mulin
		(Element& r, const Element& a) const
	{
		r = Element(Compute_t(r)*Compute_t(a) % Compute_t(_p));
		return r;
	}

	template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::divin
		(Element& r, const Element& a) const
	{
		Element ia;
		return mulin(r, inv(ia, a));
	}

	template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::addin
		(Element& r, const Element& a) const
	{
		int32_t tmp = r;
		__GIVARO_MODULAR_INTEGER_ADDIN(tmp,_p, a);
		return r = Element(tmp);
	}

	template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::subin
		(Element& r, const Element& a) const
	{
		int32_t tmp = r;
		__GIVARO_MODULAR_INTEGER_SUBIN(tmp,_p, a);
		return r = (Modular<int32_t, COMP>::Element)tmp;
	}

	template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::negin
		(Element& r) const
	{
		return __GIVARO_MODULAR_INTEGER_NEGIN(r,_p);
	}

	template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::invin
		(Element& r) const
	{
		return inv(r, r);
	}
	
	template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::axpy
		(Element& r, const Element& a, const Element& b, const Element& c) const
	{
		return __GIVARO_MODULAR_INTEGER_MULADD(r, _p, a, b, c);
	}

	template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element&  Modular<int32_t, COMP>::axpyin
		(Element& r, const Element& a, const Element& b) const
	{
		return __GIVARO_MODULAR_INTEGER_MULADDIN(r, _p, a, b);
	}
	
	template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::maxpy
		(Element& r, const Element& a, const Element& b, const Element& c) const
	{
		int32_t tmp;
		__GIVARO_MODULAR_INTEGER_MUL(tmp, _p, a, b);
		__GIVARO_MODULAR_INTEGER_SUB(r, _p, c, tmp);
		return r;
	}

	template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element&  Modular<int32_t, COMP>::axmy
		(Element& r, const Element& a, const Element& b, const Element& c) const
	{
		int32_t tmp;
		__GIVARO_MODULAR_INTEGER_MULSUB(tmp, _p, a, b, c);
		return r = (Modular<int32_t, COMP>::Element)tmp;
	}

	template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element&  Modular<int32_t, COMP>::maxpyin
		(Element& r, const Element& a, const Element& b) const
	{
		__GIVARO_MODULAR_INTEGER_SUBMULIN(r, _p, a, b );
		return r;
	}

	template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element&  Modular<int32_t, COMP>::axmyin
		(Element& r, const Element& a, const Element& b) const
	{
		return __GIVARO_MODULAR_INTEGER_MULSUB(r, _p, a, b, r );
	}
	
	// ----------------------------------
	// ----- Classic arithmetic on arrays

	template<typename COMP>
    inline void Modular<int32_t, COMP>::mul
		(const size_t sz, Array r, constArray a, constArray b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int32_t tmp;
			__GIVARO_MODULAR_INTEGER_MUL(tmp, _p,a[i], b[i]);
			r[i] = (Modular<int32_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int32_t, COMP>::mul
		(const size_t sz, Array r, constArray a, Element b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int32_t tmp;
			__GIVARO_MODULAR_INTEGER_MUL(tmp, _p, a[i], b);
			r[i] = (Modular<int32_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int32_t, COMP>::div
		(const size_t sz, Array r, constArray a, constArray b) const
	{
		for ( size_t i=sz ; --i ; ) {
			div( r[i], a[i], b[i]);
		}
	}

	template<typename COMP>
    inline void Modular<int32_t, COMP>::div
		(const size_t sz, Array r, constArray a, Element b) const
	{
		Modular<int32_t, COMP>::Element ib;
		inv(ib, b);
		mul(sz, r, a, ib);
	}

	template<typename COMP>
    inline void Modular<int32_t, COMP>::add
		(const size_t sz, Array r, constArray a, constArray b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int32_t tmp;
			__GIVARO_MODULAR_INTEGER_ADD(tmp, _p, a[i], b[i]);
			r[i] = (Modular<int32_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int32_t, COMP>::add
		(const size_t sz, Array r, constArray a, Element b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int32_t tmp;
			__GIVARO_MODULAR_INTEGER_ADD(tmp,_p, a[i], b);
			r[i] = (Modular<int32_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int32_t, COMP>::sub
		(const size_t sz, Array r, constArray a, constArray b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int32_t tmp;
			__GIVARO_MODULAR_INTEGER_SUB(tmp, _p, a[i], b[i]);
			r[i] = (Modular<int32_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int32_t, COMP>::sub
		(const size_t sz, Array r, constArray a, Element b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int32_t tmp;
			__GIVARO_MODULAR_INTEGER_SUB(tmp, _p, a[i], b);
			r[i] = (Modular<int32_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int32_t, COMP>::neg
		(const size_t sz, Array r, constArray a) const
	{
		for ( size_t i=sz ; --i ; ) {
			int32_t tmp;
			__GIVARO_MODULAR_INTEGER_NEG(tmp, _p, a[i]);
			r[i] = (Modular<int32_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int32_t, COMP>::axpy
		(const size_t sz, Array r, constArray a, constArray x, constArray y) const
	{
		for ( size_t i=sz ; --i ; ) {
			int32_t tmp;
			__GIVARO_MODULAR_INTEGER_MULADD(tmp, _p, a[i], x[i], y[i]);
			r[i] = (Modular<int32_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int32_t, COMP>::axpyin
		(const size_t sz, Array r, constArray a, constArray x) const
	{
		for ( size_t i=sz ; --i ; ) {
			int32_t tmp = r[i];
			__GIVARO_MODULAR_INTEGER_MULADDIN(tmp, _p, a[i], x[i]);
			r[i] = (Modular<int32_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int32_t, COMP>::axmy
		(const size_t sz, Array r, constArray a, constArray x, constArray y) const
	{
		for ( size_t i=sz; i--; ) {
			int32_t tmp;
			__GIVARO_MODULAR_INTEGER_MULSUB(tmp, _p, a[i], x[i], y[i]);
			r[i] = (Modular<int32_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int32_t, COMP>::maxpyin
		(const size_t sz, Array r, constArray a, constArray x) const
	{
		for ( size_t i=sz ; --i ; ) {
			int32_t tmp = r[i];
			__GIVARO_MODULAR_INTEGER_SUBMULIN(tmp, _p, a[i], x[i]);
			r[i] = (Modular<int32_t, COMP>::Element)tmp;
		}
	}
	
	// --------------------
	// ----- Initialisation

	template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::init ( Element& r, const double a ) const
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
		if (r && (sign ==-1)) r = int32_t(_p) - r;
		return r;
	}

	template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element&  Modular<int32_t, COMP>::init ( Element& r, const float a ) const
	{
		return init(r, (double)a);
	}



	template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element&  Modular<int32_t, COMP>::init ( Element& r, const unsigned long a ) const
	{
	       	return r = (Element)( a >= (unsigned long)_p ? a % (unsigned long)_p : a);
	}

	template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element&  Modular<int32_t, COMP>::init ( Element& r, const long a ) const
	{
		int sign;
		uint32_t ua;
		if (a <0) {
			sign =-1;
			ua = uint32_t(-a);
		}
		else {
			ua = uint32_t(a);
		       	sign = 1;
		}
		r = Element((ua >=_p)? ua % _p : ua);
		if (r && (sign == -1))
			r = int32_t(_p) - r;
		return r;
	}

	template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element&  Modular<int32_t, COMP>::init ( Element& r, const Integer& Residu ) const
	{
		int32_t tr;
		if (Residu <0) {
			// -a = b [p]
			// a = p-b [p]
			if ( Residu <= (Integer)(-_p) )
			       	tr = int32_t( (-Residu) % _p) ;
			else
				tr = int32_t(-Residu);
			if (tr)
				return r = Element(_p - uint32_t(tr));
			else
				return r = zero;
		}
		else {
			if (Residu >= (Integer)_p ) tr = int32_t(Residu % _p) ;
			else tr = int32_t(Residu);
			return r = Element(tr);
		}
	}

	template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::init( Element& a, const int i) const
	{
		return init(a,(long)i);
	}

	template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::init( Element& a, const unsigned int i) const
	{
		return init(a,(unsigned long)i);
	}


	template<typename COMP>
    inline void Modular<int32_t, COMP>::assign
	( const size_t sz, Array r, constArray a ) const
	{
		for ( size_t i=sz ; --i ; )
			r[i] = a[i];
	}

	template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element&  Modular<int32_t, COMP>::assign ( Element& r, const Element a ) const
	{
		return r=a;
	}

	template<typename COMP>
    inline void Modular<int32_t, COMP>::init
	( const size_t sz, Array r, constArray a ) const
	{
		for ( size_t i=sz ; --i ; )
			r[i] = a[i];
	}

	template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element& Modular<int32_t, COMP>::init ( Element& r ) const
	{
		return r = zero;
	}

	template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element&  Modular<int32_t, COMP>::dotprod
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

	template<typename COMP>
    inline typename Modular<int32_t, COMP>::Element&  Modular<int32_t, COMP>::dotprod
	( Element& r, const size_t sz, constArray a, constArray b ) const
	{
		return Modular<int32_t, COMP>::dotprod(r, _p, sz, a, b);
	}


	//  a -> r: int32_t to double
	template<typename COMP>
    inline void
	Modular<int32_t, COMP>::i2d ( const size_t sz, double* r, constArray a ) const
	{
		for (size_t i=0; i<sz; ++i) r[i] = a[i];
	}

	//  a -> r: double to int32_t
	template<typename COMP>
    inline void
	Modular<int32_t, COMP>::d2i ( const size_t sz, Array r, const double* a ) const
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
			if (r[i] <_p) r[i] %= _p;
		}
		//    r[i] = (tmp.r[1] <_p ? tmp.r[1] : tmp.r[1]-_p);
		//    r[i] = (r[i] <_p ? r[i] : r[i]%_p);
		//    r[i] = (tmp.r[1] <_p ? tmp.r[1] : tmp.r[1]%_p);
	}



	// -- Input: (z, <_p>)
	template<typename COMP>
    inline std::istream& Modular<int32_t, COMP>::read (std::istream& s)
	{
		char ch;
		s >> std::ws >> ch;
		if (ch != '(')
			//    GivError::throw_error( GivBadFormat("Modular<int32_t, COMP>::read: syntax error: no '('"));
			std::cerr << "GivBadFormat(Modular<int32_t, COMP>::read: syntax error: no '('))" << std::endl;

		s >> std::ws >> ch;
		if (ch != 'z')
			//    GivError::throw_error( GivBadFormat("Modular<int32_t, COMP>::read: bad domain object"));
			std::cerr << "GivBadFormat(Modular<int32_t, COMP>::read: bad domain object))" << std::endl;

		s >> std::ws >> ch;
		if (ch != ',')
			//    GivError::throw_error( GivBadFormat("Modular<int32_t, COMP>::read: syntax error: no ','"));
			std::cerr << "GivBadFormat(Modular<int32_t, COMP>::read: syntax error: no ',')) " << std::endl;

		s >> std::ws >> _p;

		s >> std::ws >> ch;
		if (ch != ')')
			//    GivError::throw_error( GivBadFormat("Modular<int32_t, COMP>::read: syntax error: no ')'"));
			std::cerr << "GivBadFormat(Modular<int32_t, COMP>::read: syntax error: no ')')) " << std::endl;

		return s;
	}

	template<>
	inline std::ostream& Modular<int32_t, int32_t>::write (std::ostream& s) const
	{
		return s << "Modular<int32_t, uint32_t> modulo " << residu();
	}

	template<>
	inline std::ostream& Modular<int32_t, uint32_t>::write (std::ostream& s) const
	{
		return s << "Modular<int32_t, uint32_t> modulo " << residu();
	}

	template<>
	inline std::ostream& Modular<int32_t, int64_t>::write (std::ostream& s) const
	{
		return s << "Modular<int32_t, uint64_t> modulo " << residu();
	}

	template<>
	inline std::ostream& Modular<int32_t, uint64_t>::write (std::ostream& s) const
	{
		return s << "Modular<int32_t, uint64_t> modulo " << residu();
	}

	template<typename COMP>
    inline std::istream& Modular<int32_t, COMP>::read (std::istream& s, Element& a) const
	{
		Integer tmp;
		s >> tmp;
		init(a, tmp);
		return s;
	}

	template<typename COMP>
	inline std::ostream& Modular<int32_t, COMP>::write (std::ostream& s, const Element a) const
	{
		return s << a;
	}

} // namespace Givaro

#endif // __GIVARO_zpz32std_INL
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
