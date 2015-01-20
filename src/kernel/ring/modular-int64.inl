// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givzpz64std.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givzpz64std.inl,v 1.19 2011-02-02 17:16:43 bboyer Exp $
// ==========================================================================

#ifndef __GIVARO_zpz64std_INL
#define __GIVARO_zpz64std_INL

#include <givaro/modular-defines.h>

namespace Givaro {

	// ------------------------
	// ----- Classic arithmetic
	
	inline Modular<int64_t>::Element& Modular<int64_t>::mul
		(Element& r, const Element& a, const Element& b) const
	{
		return  __GIVARO_MODULAR_INTEGER_MUL(r,(int64_t)_p,(int64_t)a,(int64_t)b);
	}

	inline Modular<int64_t>::Element& Modular<int64_t>::sub
		(Element& r, const Element& a, const Element& b) const
	{
		return __GIVARO_MODULAR_INTEGER_SUB(r,(int64_t)_p,(int64_t)a,(int64_t)b);
	}

	inline Modular<int64_t>::Element& Modular<int64_t>::add
		(Element& r, const Element& a, const Element& b) const
	{
		__GIVARO_MODULAR_INTEGER_ADD(r,(int64_t)_p,(int64_t)a,(int64_t)b);
		return r;
	}

	inline Modular<int64_t>::Element& Modular<int64_t>::neg
		(Element& r, const Element& a) const
	{
		return __GIVARO_MODULAR_INTEGER_NEG(r,(int64_t)_p,(int64_t)a);
	}

	inline Modular<int64_t>::Element& Modular<int64_t>::inv
		(Element& r, const Element& a) const
	{
		int64_t u;
		Modular<int64_t>::invext(u, a, (int64_t)_p);
		return r = (u<0)?(Modular<int64_t>::Element)u + (int64_t)_p:(Modular<int64_t>::Element)u;
	}

	inline Modular<int64_t>::Element& Modular<int64_t>::div
		(Element& r, const Element& a, const Element& b) const
	{
		return mulin( inv(r,b), a );
	}
	inline Modular<int64_t>::Element& Modular<int64_t>::mulin
		(Element& r, const Element& a) const
	{
		return __GIVARO_MODULAR_INTEGER_MULIN(r,(int64_t)_p, (int64_t)a);
	}

	inline Modular<int64_t>::Element& Modular<int64_t>::divin
		(Element& r, const Element& a) const
	{
		Modular<int64_t>::Element ia;
		inv(ia, a);
		return mulin(r, ia);
	}

	inline Modular<int64_t>::Element& Modular<int64_t>::addin
		(Element& r, const Element& a) const
	{
		int64_t tmp = (int64_t)r;
		__GIVARO_MODULAR_INTEGER_ADDIN(tmp,(int64_t)_p, (int64_t)a);
		return r = (Modular<int64_t>::Element)tmp;
	}

	inline Modular<int64_t>::Element& Modular<int64_t>::subin
		(Element& r, const Element& a) const
	{
		int64_t tmp = (int64_t)r;
		__GIVARO_MODULAR_INTEGER_SUBIN(tmp,(int64_t)_p, (int64_t)a);
		return r = (Modular<int64_t>::Element)tmp;
	}

	inline Modular<int64_t>::Element& Modular<int64_t>::negin
		(Element& r) const
	{
		return __GIVARO_MODULAR_INTEGER_NEGIN(r,(int64_t)_p);
	}

	inline Modular<int64_t>::Element& Modular<int64_t>::invin
		(Element& r) const
	{
		int64_t u;
		Modular<int64_t>::invext(u, r, (int64_t)_p);
		return r = (u<0)?(Modular<int64_t>::Element)u + (int64_t)_p:(Modular<int64_t>::Element)u;
	}
	
	inline Modular<int64_t>::Element& Modular<int64_t>::axpy
		(Element& r, const Element& a, const Element& b, const Element& c) const
	{
		return __GIVARO_MODULAR_INTEGER_MULADD(r, (int64_t)_p, (int64_t)a, (int64_t)b, (int64_t)c);
	}

	inline Modular<int64_t>::Element&  Modular<int64_t>::axpyin
		(Element& r, const Element& a, const Element& b) const
	{
		return __GIVARO_MODULAR_INTEGER_MULADDIN(r, (int64_t)_p, (int64_t)a, (int64_t)b);
	}
	
	inline Modular<int64_t>::Element& Modular<int64_t>::maxpy
		(Element& r, const Element& a, const Element& b, const Element& c) const
	{
		int64_t tmp;
		__GIVARO_MODULAR_INTEGER_MUL(tmp, (int64_t)_p, (int64_t)a, (int64_t)b);
		__GIVARO_MODULAR_INTEGER_SUB(r, (int64_t)_p, (int64_t)c, tmp);
		return r;
	}

	inline Modular<int64_t>::Element&  Modular<int64_t>::axmy
		(Element& r, const Element& a, const Element& b, const Element& c) const
	{
		int64_t tmp;
		__GIVARO_MODULAR_INTEGER_MULSUB(tmp, (int64_t)_p, (int64_t)a, (int64_t)b, (int64_t)c);
		return r = (Modular<int64_t>::Element)tmp;
	}

	inline Modular<int64_t>::Element&  Modular<int64_t>::maxpyin
		(Element& r, const Element& a, const Element& b) const
	{
		__GIVARO_MODULAR_INTEGER_SUBMULIN(r, (int64_t)_p, (int64_t)a, (int64_t)b );
		return r;
	}

	inline Modular<int64_t>::Element&  Modular<int64_t>::axmyin
		(Element& r, const Element& a, const Element& b) const
	{
		return __GIVARO_MODULAR_INTEGER_MULSUB(r, (int64_t)_p, (int64_t)a, (int64_t)b, r );
	}
	
	// ----------------------------------
	// ----- Classic arithmetic on arrays

	inline void Modular<int64_t>::mul
		(const size_t sz, Array r, constArray a, constArray b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int64_t tmp;
			__GIVARO_MODULAR_INTEGER_MUL(tmp, (int64_t)_p,(int64_t)a[i], (int64_t)b[i]);
			r[i] = (Modular<int64_t>::Element)tmp;
		}
	}

	inline void Modular<int64_t>::mul
		(const size_t sz, Array r, constArray a, Element b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int64_t tmp;
			__GIVARO_MODULAR_INTEGER_MUL(tmp, (int64_t)_p, (int64_t)a[i], (int64_t)b);
			r[i] = (Modular<int64_t>::Element)tmp;
		}
	}

	inline void Modular<int64_t>::div
		(const size_t sz, Array r, constArray a, constArray b) const
	{
		for ( size_t i=sz ; --i ; ) {
			div( r[i], a[i], b[i]);
		}
	}

	inline void Modular<int64_t>::div
		(const size_t sz, Array r, constArray a, Element b) const
	{
		Modular<int64_t>::Element ib;
		inv(ib, b);
		mul(sz, r, a, ib);
	}

	inline void Modular<int64_t>::add
		(const size_t sz, Array r, constArray a, constArray b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int64_t tmp;
			__GIVARO_MODULAR_INTEGER_ADD(tmp, (int64_t)_p, (int64_t)a[i], (int64_t)b[i]);
			r[i] = (Modular<int64_t>::Element)tmp;
		}
	}

	inline void Modular<int64_t>::add
		(const size_t sz, Array r, constArray a, Element b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int64_t tmp;
			__GIVARO_MODULAR_INTEGER_ADD(tmp,(int64_t)_p, (int64_t)a[i], (int64_t)b);
			r[i] = (Modular<int64_t>::Element)tmp;
		}
	}

	inline void Modular<int64_t>::sub
		(const size_t sz, Array r, constArray a, constArray b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int64_t tmp;
			__GIVARO_MODULAR_INTEGER_SUB(tmp, (int64_t)_p, (int64_t)a[i], (int64_t)b[i]);
			r[i] = (Modular<int64_t>::Element)tmp;
		}
	}

	inline void Modular<int64_t>::sub
		(const size_t sz, Array r, constArray a, Element b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int64_t tmp;
			__GIVARO_MODULAR_INTEGER_SUB(tmp, (int64_t)_p, (int64_t)a[i], (int64_t)b);
			r[i] = (Modular<int64_t>::Element)tmp;
		}
	}

	inline void Modular<int64_t>::neg
		(const size_t sz, Array r, constArray a) const
	{
		for ( size_t i=sz ; --i ; ) {
			int64_t tmp;
			__GIVARO_MODULAR_INTEGER_NEG(tmp, (int64_t)_p, (int64_t)a[i]);
			r[i] = (Modular<int64_t>::Element)tmp;
		}
	}

	inline void Modular<int64_t>::axpy
		(const size_t sz, Array r, constArray a, constArray x, constArray y) const
	{
		for ( size_t i=sz ; --i ; ) {
			int64_t tmp;
			__GIVARO_MODULAR_INTEGER_MULADD(tmp, (int64_t)_p, (int64_t)a[i], (int64_t)x[i], (int64_t)y[i]);
			r[i] = (Modular<int64_t>::Element)tmp;
		}
	}

	inline void Modular<int64_t>::axpyin
		(const size_t sz, Array r, constArray a, constArray x) const
	{
		for ( size_t i=sz ; --i ; ) {
			int64_t tmp = (int64_t)r[i];
			__GIVARO_MODULAR_INTEGER_MULADDIN(tmp, (int64_t)_p, (int64_t)a[i], (int64_t)x[i]);
			r[i] = (Modular<int64_t>::Element)tmp;
		}
	}

	inline void Modular<int64_t>::axmy
		(const size_t sz, Array r, constArray a, constArray x, constArray y) const
	{
		for ( size_t i=sz; i--; ) {
			int64_t tmp;
			__GIVARO_MODULAR_INTEGER_MULSUB(tmp, (int64_t)_p, (int64_t)a[i], (int64_t)x[i], (int64_t)y[i]);
			r[i] = (Modular<int64_t>::Element)tmp;
		}
	}

	inline void Modular<int64_t>::maxpyin
		(const size_t sz, Array r, constArray a, constArray x) const
	{
		for ( size_t i=sz ; --i ; ) {
			int64_t tmp = (int64_t)r[i];
			__GIVARO_MODULAR_INTEGER_SUBMULIN(tmp, (int64_t)_p, (int64_t)a[i], (int64_t)x[i]);
			r[i] = (Modular<int64_t>::Element)tmp;
		}
	}
	
	// --------------------
	// ----- Initialisation
	
	inline  Modular<int64_t>::Element&  Modular<int64_t>::init ( Element& r, const unsigned long a ) const
	{
		return r = (Element)( a >= (uint64_t)_p ? a % (uint64_t)_p : a);
	}

	inline  Modular<int64_t>::Element&  Modular<int64_t>::init ( Element& r, const long a ) const
	{
		int64_t sign; uint64_t ua;
		if (a <0) { sign =-1; ua = (unsigned int)-a;}
		else { ua = (unsigned int)a; sign =1; }
		r = (Element)((ua >=_p) ? ua % (uint64_t)_p : ua);
		if (r && (sign ==-1)) r = (Element)_p - r;
		return r;
	}


	inline Modular<int64_t>::Element&  Modular<int64_t>::init ( Element& r, const Integer& Residu ) const
	{
		if (Residu <0) {
			int64_t tr;
			// -a = b [p]
			// a = p-b [p]
			if ( (-Residu) >= (Integer)(_p) ) tr = int64_t( (-Residu) % (Integer)_p) ;
			else tr = int64_t(-Residu);
			if (tr) return r = (Element)( (uint64_t)_p - (uint64_t)tr ) ;
			else return r = zero;
		} else {
			Integer ip(_p);
			if (Residu >= ip ) return r =   int64_t(Residu % ip) ;
			else return r = int64_t(Residu);
		}
	}

	inline  Modular<int64_t>::Element& Modular<int64_t>::init( Element& a, const int i) const
	{
		return init(a,(long)i);
	}
	
	inline  Modular<int64_t>::Element& Modular<int64_t>::init( Element& a, const unsigned int i) const
	{
		return init(a,(unsigned long)i);
	}


	inline  Modular<int64_t>::Element&  Modular<int64_t>::init ( Element& r, const unsigned long long a ) const
	{
		return r = (Element)( a >= (uint64_t)_p ? a % (uint64_t)_p : a);
	}

	inline  Modular<int64_t>::Element&  Modular<int64_t>::init ( Element& r, const double a ) const
	{
		return init(r, (int64_t)a);
	}

	inline  Modular<int64_t>::Element&  Modular<int64_t>::init ( Element& r, const float a ) const
	{
		return init(r, (double)a);
	}

	inline  Modular<int64_t>::Element&  Modular<int64_t>::init ( Element& r, const long long a ) const
	{
		int sign; uint64_t ua;
		if (a <0) { sign =-1; ua = (unsigned int)-a;}
		else { ua = (unsigned int)a; sign =1; }
		r = (Element) ( (ua >=_p) ? ua % (uint64_t)_p : ua) ;
		if (r && (sign ==-1)) r = (Element)_p - r;
		return r;
	}

	inline void Modular<int64_t>::init ( const size_t sz, Array r, constArray a ) const
	{
		for ( size_t i = sz; --i; ) {
			r[i] = a[i];
		}
	}

	inline Modular<int64_t>::Element& Modular<int64_t>::init ( Element& r ) const
	{
		return r = zero;
	}
	
	// ---------
	// -- Assign

	inline void Modular<int64_t>::assign( const size_t sz, Array r, constArray a ) const
	{
		for ( size_t i = sz; --i; )
		    r[i] = a[i];
	}

	inline Modular<int64_t>::Element&  Modular<int64_t>::assign
	  ( Element& r, const Element a ) const
	{
		return r = a;
	}
	
	// ----------
	// -- Convert
	
	template<class XXX> inline XXX& Modular<int64_t>::convert(XXX& s, const Element a) const
	{
		return s = XXX(a);
	}
	
	template<> inline Integer& Modular<int64_t>::convert(Integer& i, const Element a) const
	{
		unsigned long ur;
		return i = (Integer)convert(ur, a);
	}






inline Modular<int64_t>::Element& Modular<int64_t>::dotprod
  ( Element& r, const int bound, const size_t sz, constArray a, constArray b ) const
{
  unsigned int stride = 1;
  if ((int64_t)bound < Signed_Trait<Element>::max() )
   stride = (unsigned int) ( GIVARO_MAXULONG/((unsigned long)bound * (unsigned long)bound) );
  unsigned long dot = (unsigned long) zero; // this is intented !
  if ((sz <10) && (sz <stride)) {
    for(  size_t i= sz; i--; )
#ifdef __x86_64__
      dot += (unsigned long)a[i] * (unsigned long)b[i];
#else
      dot = (unsigned long) (dot + a[i] * b[i]);
#endif
    if (dot > _p) r = (Element)(dot % (uint64_t)_p);
    else r = (Element)dot;
    return r;
  }
  unsigned int i_begin=0;
  stride &= (unsigned int)~0x1;
  if (stride ==0) {
    for(  size_t i= sz; --i; ) {
#ifdef __x86_64__
      dot += (unsigned long)a[i] * (unsigned long)b[i];
      if (dot>_p) dot %= _p;
#else
      dot = (unsigned long) (dot + a[i] * b[i]);
      if (dot>_p) dot = (unsigned long) (dot % _p);
#endif
    }
    r = (Element)dot;
    return r;
  }
  do {
    size_t min_sz = ((sz-i_begin) < stride ? (sz-i_begin) : stride);
    if ((min_sz & 0x1) !=0)
      {
	      --min_sz;
	      ++i_begin;
#ifdef __x86_64__
	      dot += (unsigned long)a++[min_sz] * (unsigned long)b++[min_sz];
#else
	      dot = (unsigned long) (dot + a++[min_sz] * b++[min_sz]);
#endif
      }
    if (min_sz > 1)
      for(  size_t i= min_sz; i>0; --i, --i, ++a, ++a, ++b, ++b ) //!@todo o_O
      {
#ifdef __x86_64__
        dot += (unsigned long)a[0] * (unsigned long)b[0];
        dot += (unsigned long)a[1] * (unsigned long)b[1];
#else
	dot = (unsigned long) (dot +  a[0] * b[0] );
	dot = (unsigned long) (dot +  a[1] * b[1] );
#endif

      }
#ifdef __x86_64__
    if (dot>(uint64_t)_p) dot %= (uint64_t)_p;
#else
    if (dot>_p) dot = (unsigned long) (dot % _p);
#endif
    i_begin += (unsigned int) min_sz;
  } while (i_begin <sz);
  r = (Element)dot;
  return r;
}

inline Modular<int64_t>::Element& Modular<int64_t>::dotprod
  ( Element& r, const size_t sz, constArray a, constArray b ) const
{
	return Modular<int64_t>::dotprod(r, int(_p), sz, a, b);
}


  //  a -> r: int64_t to double
inline void
  Modular<int64_t>::i2d ( const size_t sz, double* r, constArray a ) const
{
  for (size_t i=0; i<sz; ++i)  {
	  r[i] = (double) a[i];
  }
}

  //  a -> r: double to int64_t
inline void
  Modular<int64_t>::d2i ( const size_t sz, Array r, const double* a ) const
{
  union d_2_l {
    double d;
    int64_t r[2];
  };
//  static const double offset = 4503599627370496.0; // 2^52
  double offset = 4503599627370496.0; // 2^52
  for (size_t i=0; i<sz; ++i)
  {
       d_2_l tmp;
      // - normalization: put fractional part at the end of the representation
      tmp.d = a[i] + offset;
      r[i] = tmp.r[1];
      if (r[i] <(int64_t)_p) r[i] %= _p;
  }
  //    r[i] = (tmp.r[1] <_p ? tmp.r[1] : tmp.r[1]-_p);
  //    r[i] = (r[i] <_p ? r[i] : r[i]%_p);
  //    r[i] = (tmp.r[1] <_p ? tmp.r[1] : tmp.r[1]%_p);
}


 // -- Input: (z, <_p>)
inline std::istream& Modular<int64_t>::read (std::istream& s)
{
  char ch;
  s >> std::ws >> ch;
  if (ch != '(')
//    GivError::throw_error( GivBadFormat("Modular<int64_t>::read: syntax error: no '('"));
    std::cerr << "GivBadFormat(Modular<int64_t>::read: syntax error: no '('))" << std::endl;

  s >> std::ws >> ch;
  if (ch != 'z')
//    GivError::throw_error( GivBadFormat("Modular<int64_t>::read: bad domain object"));
    std::cerr << "GivBadFormat(Modular<int64_t>::read: bad domain object))" << std::endl;

  s >> std::ws >> ch;
  if (ch != ',')
//    GivError::throw_error( GivBadFormat("Modular<int64_t>::read: syntax error: no ','"));
    std::cerr << "GivBadFormat(Modular<int64_t>::read: syntax error: no ',')) " << std::endl;

  s >> std::ws >> _p;


  s >> std::ws >> ch;
  if (ch != ')')
//    GivError::throw_error( GivBadFormat("Modular<int64_t>::read: syntax error: no ')'"));
    std::cerr << "GivBadFormat(Modular<int64_t>::read: syntax error: no ')')) " << std::endl;

  return s;
}

inline std::ostream& Modular<int64_t>::write (std::ostream& s ) const
{
  	return s << "Modular<int64_t> modulo " << residu();
}

inline std::istream& Modular<int64_t>::read (std::istream& s, Element& a) const
{
	Integer tmp;
	s >> tmp;
	init(a, tmp);
	return s;
}

inline std::ostream& Modular<int64_t>::write (std::ostream& s, const Element a) const
{
  	return s << a;
}


} // namespace Givaro

#endif // __GIVARO_zpz64std_INL
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
