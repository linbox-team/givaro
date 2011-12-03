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

// Description:

// ---------
// -- normalized operations
// ---------

#include <givaro/givconfig.h>

#ifndef __DONOTUSE_Givaro_SIXTYFOUR__

// r = a*b
#define __GIVARO_ZPZ64_N_MUL(r,p,a,b) ( r = (Rep) ( (uint64_t)(a*b) % (uint64_t)p ) )
// r *= a
#define __GIVARO_ZPZ64_N_MULIN(r,p,a) (  r = (Rep) ( (uint64_t)(r*a) % (uint64_t)p  ) )

// r = a - b
#define __GIVARO_ZPZ64_N_SUB(r,p,a,b) { r = (a-b); r= (r < 0 ? r+p : r);}
// r -= a
#define __GIVARO_ZPZ64_N_SUBIN(r,p,a) { r -= a; r= (r < 0 ? r+p : r);}

// r = a+b
#define __GIVARO_ZPZ64_N_ADD(r,p,a,b) { r = (a+b); r= (r < p ? r : r-p);}
// r += a
#define __GIVARO_ZPZ64_N_ADDIN(r,p,a) { r += a;  r= (r < p ? r : r-p);}

// r <- a*b+c % p
#define __GIVARO_ZPZ64_N_MULADD(r,p,a,b,c) ( r = (Rep) ( (uint64_t)(a*b+c) % (uint64_t)p ) )

#define __GIVARO_ZPZ64_N_MULADDIN(r,p,a,b) ( r = (Rep) ( (uint64_t)(a*b+r) % (uint64_t)p ) )

#define __GIVARO_ZPZ64_N_NEG(r,p,a) { r = (a == 0 ? 0 : p-a); }
#define __GIVARO_ZPZ64_N_NEGIN(r,p) { r = (r == 0 ? 0 : p-r); }

// a*b-c
#define __GIVARO_ZPZ64_N_MULSUB(r,p,a,b,c) ( r = (Rep) ( (uint64_t)(a*b+p-c) % (uint64_t)p ) )

// r-a*b
#define __GIVARO_ZPZ64_N_SUBMULIN(r,p,a,b) { \
    __GIVARO_ZPZ64_N_MULSUB(r,p,a,b,r); __GIVARO_ZPZ64_N_NEGIN(r,p); }

namespace Givaro {

inline ZpzDom<Std64>::Residu_t ZpzDom<Std64>::residu( ) const
{ return _p; }

inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::mul (Rep& r, const Rep a, const Rep b) const
{
   int64_t tmp;
  __GIVARO_ZPZ64_N_MUL(tmp,(int64_t)_p,(int64_t)a,(int64_t)b);
  return r = (ZpzDom<Std64>::Rep)tmp;
}

inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::neg (Rep& r, const Rep a) const
{
   int64_t tmp;
  __GIVARO_ZPZ64_N_NEG(tmp,(int64_t)_p,(int64_t)a);
  return r = (ZpzDom<Std64>::Rep)tmp;
}

inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::inv (Rep& r, const Rep a) const
{
//    int64_t d, u, v;
   int64_t u;
  ZpzDom<Std64>::invext(u, a, (int64_t)_p);
  return r = (u<0)?(ZpzDom<Std64>::Rep)(u+(int64_t)_p):(ZpzDom<Std64>::Rep)u;
 }

inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::div (Rep& r, const Rep a, const Rep b) const
{
   int64_t tmp;
   int64_t ib;
  inv(ib, b);
  __GIVARO_ZPZ64_N_MUL(tmp,(int64_t)_p,(int64_t)a,(int64_t)ib);
  return r = (ZpzDom<Std64>::Rep)tmp;
}

inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::sub (Rep& r, const Rep a, const Rep b) const
{
   int64_t tmp;
  __GIVARO_ZPZ64_N_SUB(tmp,(int64_t)_p,(int64_t)a,(int64_t)b);
  return r = (ZpzDom<Std64>::Rep)tmp;
}

inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::add (Rep& r, const Rep a, const Rep b) const
{
   int64_t tmp;
  __GIVARO_ZPZ64_N_ADD(tmp,(int64_t)_p,(int64_t)a,(int64_t)b);
  return r = (ZpzDom<Std64>::Rep)tmp;
}

 // -- inline array operations between ZpzDom<Std64>::Rep
inline void ZpzDom<Std64>::mul (const size_t sz, Array r, constArray a, constArray b) const
{
	for ( size_t i = sz; --i; ) {
		int64_t tmp;
		__GIVARO_ZPZ64_N_MUL(tmp, (int64_t)_p,(int64_t)a[i], (int64_t)b[i]);
		r[i] = (ZpzDom<Std64>::Rep)tmp;
	}
}

inline void ZpzDom<Std64>::mul (const size_t sz, Array r, constArray a, Rep b) const
{
	for ( size_t i = sz; --i; ) {
		int64_t tmp;
		__GIVARO_ZPZ64_N_MUL(tmp, (int64_t)_p, (int64_t)a[i], (int64_t)b);
		r[i] = (ZpzDom<Std64>::Rep)tmp;
	}
}

inline void ZpzDom<Std64>::div (const size_t sz, Array r, constArray a, constArray b) const
{
	for ( size_t i = sz; --i; ) {
		div( r[i], a[i], b[i]);
	}
}

inline void ZpzDom<Std64>::div (const size_t sz, Array r, constArray a, Rep b) const
{
  ZpzDom<Std64>::Rep ib;
  inv(ib, b);
  mul(sz, r, a, ib);
}

inline void ZpzDom<Std64>::add (const size_t sz, Array r, constArray a, constArray b) const
{
	for ( size_t i = sz; --i; ) {
		int64_t tmp;
		__GIVARO_ZPZ64_N_ADD(tmp, (int64_t)_p, (int64_t)a[i], (int64_t)b[i]);
		r[i] = (ZpzDom<Std64>::Rep)tmp;
	}
}

inline void ZpzDom<Std64>::add (const size_t sz, Array r, constArray a, Rep b) const
{
	for ( size_t i = sz; --i; ) {
		int64_t tmp;
		__GIVARO_ZPZ64_N_ADD(tmp, (int64_t)_p, (int64_t)a[i], (int64_t)b);
		r[i] = (ZpzDom<Std64>::Rep)tmp;
	}
}

inline void ZpzDom<Std64>::sub (const size_t sz, Array r, constArray a, constArray b) const
{
	for ( size_t i = sz; --i; ) {
		int64_t tmp;
		__GIVARO_ZPZ64_N_SUB(tmp, (int64_t)_p, (int64_t)a[i], (int64_t)b[i]);
		r[i] = (ZpzDom<Std64>::Rep)tmp;
	}
}

inline void ZpzDom<Std64>::sub (const size_t sz, Array r, constArray a, Rep b) const
{
	for ( size_t i = sz; --i; ) {
		int64_t tmp;
		__GIVARO_ZPZ64_N_SUB(tmp, (int64_t)_p, (int64_t)a[i], (int64_t)b);
		r[i] = (ZpzDom<Std64>::Rep)tmp;
	}
}

inline void ZpzDom<Std64>::neg (const size_t sz, Array r, constArray a) const
{
	for ( size_t i = sz; --i; ) {
		int64_t tmp;
		__GIVARO_ZPZ64_N_NEG(tmp, (int64_t)_p, (int64_t)a[i]);
		r[i] = (ZpzDom<Std64>::Rep)tmp;
	}
}


inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::mulin (Rep& r, const Rep a) const
{
   int64_t tmp = (int64_t)r;
  __GIVARO_ZPZ64_N_MULIN(tmp,(int64_t)_p, (int64_t)a);
  return r = (ZpzDom<Std64>::Rep)tmp;
}

inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::divin (Rep& r, const Rep a) const
{
  ZpzDom<Std64>::Rep ia;
  inv(ia, a);
  return mulin(r, ia);
}

inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::addin (Rep& r, const Rep a) const
{
   int64_t tmp = (int64_t)r;
  __GIVARO_ZPZ64_N_ADDIN(tmp,(int64_t)_p, (int64_t)a);
  return r = (ZpzDom<Std64>::Rep)tmp;
}

inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::subin (Rep& r, const Rep a) const
{
   int64_t tmp = (int64_t)r;
  __GIVARO_ZPZ64_N_SUBIN(tmp,(int64_t)_p, (int64_t)a);
  return r = (ZpzDom<Std64>::Rep)tmp;
}


inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::negin (Rep& r) const
{
  __GIVARO_ZPZ64_N_NEGIN(r,(int64_t)_p);
  return r;
}

inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::invin (Rep& r) const
{
   int64_t u;
  ZpzDom<Std64>::invext(u, r, (int64_t)_p);
  return r = (u<0)?(ZpzDom<Std64>::Rep)(u+(int64_t)_p):(ZpzDom<Std64>::Rep)u;
}


inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::axpy
 (Rep& r, const Rep a, const Rep b, const Rep c) const
{
   int64_t tmp;
  __GIVARO_ZPZ64_N_MULADD(tmp, (int64_t)_p, (int64_t)a, (int64_t)b, (int64_t)c);
  return r = (ZpzDom<Std64>::Rep)tmp;
}

inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::axpyin
 (Rep& r, const Rep a, const Rep b) const
{
   int64_t tmp = (int64_t)r;
  __GIVARO_ZPZ64_N_MULADDIN(tmp, (int64_t)_p, (int64_t)a, (int64_t)b);
  return r = (ZpzDom<Std64>::Rep)tmp;
}


inline void ZpzDom<Std64>::axpy
  (const size_t sz, Array r, constArray a, constArray x, constArray y) const
{
	for ( size_t i = sz; --i; ) {
		int64_t tmp;
		__GIVARO_ZPZ64_N_MULADD(tmp, (int64_t)_p, (int64_t)a[i], (int64_t)x[i], (int64_t)y[i]);
		r[i] = (ZpzDom<Std64>::Rep)tmp;
	}
}

inline void ZpzDom<Std64>::axpyin
  (const size_t sz, Array r, constArray a, constArray x) const
{
	for ( size_t i = sz; --i; ) {
		int64_t tmp = (int64_t)r[i];
		__GIVARO_ZPZ64_N_MULADDIN(tmp, (int64_t)_p, (int64_t)a[i], (int64_t)x[i]);
		r[i] = (ZpzDom<Std64>::Rep)tmp;
	}
}

inline ZpzDom<Std64>::Rep&  ZpzDom<Std64>::axmy
 (Rep& r, const Rep a, const Rep b, const Rep c) const
{
   int64_t tmp;
  __GIVARO_ZPZ64_N_MULSUB(tmp, (int64_t)_p, (int64_t)a, (int64_t)b, (int64_t)c);
  return r = (ZpzDom<Std64>::Rep)tmp;
}

// r = c - a*b
inline ZpzDom<Std64>::Rep&  ZpzDom<Std64>::maxpy
 (Rep& r, const Rep a, const Rep b, const Rep c) const
{
   int64_t tmp = (int64_t)c;
  __GIVARO_ZPZ64_N_SUBMULIN(tmp, (int64_t)_p, (int64_t)a, (int64_t)b );
  return r = (ZpzDom<Std64>::Rep)tmp;
}

#if 0
inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::maxpy (Rep& r,
						 const Rep a, const Rep b, const Rep c) const
{
	int64_t tmp;
	__GIVARO_ZPZ64_N_MUL(tmp, (int64_t)_p, (int64_t)a, (int64_t)b);
	__GIVARO_ZPZ64_N_SUB(r, (int64_t)_p, (int64_t)c, tmp);
	return r;
}
#endif

// r -= a*b
inline ZpzDom<Std64>::Rep&  ZpzDom<Std64>::maxpyin (Rep& r, const Rep a, const Rep b) const
{
  __GIVARO_ZPZ64_N_SUBMULIN(r, (int64_t)_p, (int64_t)a, (int64_t)b );
  return r;
//    int64_t tmp = (int64_t)r;
//   __GIVARO_ZPZ64_N_SUBMULIN(tmp, (int64_t)_p, (int64_t)a, (int64_t)b );
//   return r = (ZpzDom<Std64>::Rep)tmp;
}

// r = a*b - r
inline ZpzDom<Std64>::Rep&  ZpzDom<Std64>::axmyin (Rep& r, const Rep a, const Rep b) const
{
    maxpyin(r,a,b);
    return negin(r);
}


inline void ZpzDom<Std64>::axmy (const size_t sz, Array r, constArray a, constArray x, constArray y) const
{
  for ( size_t i = sz; --i; ) {
     int64_t tmp;
    __GIVARO_ZPZ64_N_MULSUB(tmp, (int64_t)_p, (int64_t)a[i], (int64_t)x[i], (int64_t)y[i]);
    r[i] = (ZpzDom<Std64>::Rep)tmp;
  }
}

// r = a*b - r
inline void ZpzDom<Std64>::maxpyin (const size_t sz, Array r, constArray a, constArray x) const
{
	for ( size_t i = sz; --i; ) {
		int64_t tmp = (int64_t)r[i];
		__GIVARO_ZPZ64_N_SUBMULIN(tmp, (int64_t)_p, (int64_t)a[i], (int64_t)x[i]);
		r[i] = (ZpzDom<Std64>::Rep)tmp;
	}
}

 // ------------------------- Miscellaneous functions

inline int ZpzDom<Std64>::areEqual(const Rep a, const Rep b) const
{ return a == b; }

inline int ZpzDom<Std64>::areNEqual(const Rep a, const Rep b) const
{ return a != b; }

inline int ZpzDom<Std64>::isZero(const Rep a) const
{ return a == ZpzDom<Std64>::zero; }

inline int ZpzDom<Std64>::isnzero(const Rep a) const
{ return a != ZpzDom<Std64>::zero; }

inline int ZpzDom<Std64>::isOne(const Rep a) const
{ return a == ZpzDom<Std64>::one; }


inline size_t ZpzDom<Std64>::length(const Rep ) const
{ return ZpzDom<Std64>::size_rep;}


// ---------
// -- misc operations
// ---------

inline void ZpzDom<Std64>::assign
  ( const size_t sz, Array r, constArray a ) const
{
    for ( size_t i = sz; --i; )
        r[i] = a[i];
}

inline  ZpzDom<Std64>::Rep&  ZpzDom<Std64>::init ( Rep& r, const unsigned long a ) const
{ return r = (Rep)( a >= (uint64_t)_p ? a % (uint64_t)_p : a);
}


inline  ZpzDom<Std64>::Rep&  ZpzDom<Std64>::init ( Rep& r, const long a ) const
{
  int64_t sign; uint64_t ua;
  if (a <0) { sign =-1; ua = (unsigned int)-a;}
  else { ua = (unsigned int)a; sign =1; }
  r = (Rep)((ua >=_p) ? ua % (uint64_t)_p : ua);
  if (r && (sign ==-1)) r = (Rep)_p - r;
  return r;
}


inline ZpzDom<Std64>::Rep&  ZpzDom<Std64>::init ( Rep& r, const Integer& Residu ) const
{
  int64_t tr;
  if (Residu <0) {
      // -a = b [p]
      // a = p-b [p]
    if ( (-Residu) >= (Integer)(_p) ) tr = int64_t( (-Residu) % (Integer)_p) ;
    else tr = int64_t(-Residu);
    if (tr)
      return r = (Rep)( (uint64_t)_p - (uint64_t)tr ) ;
    else
      return r = zero;
  } else {
      if (Residu >= (Integer)_p ) tr =   int64_t(Residu % (Integer)_p) ;
    else tr = int64_t(Residu);
    return r = tr;
  }
}




inline  ZpzDom<Std64>::Rep& ZpzDom<Std64>::init( Rep& a, const int i) const { return init(a,(long)i); }
inline  ZpzDom<Std64>::Rep& ZpzDom<Std64>::init( Rep& a, const unsigned int i) const { return init(a,(unsigned long)i); }


inline  ZpzDom<Std64>::Rep&  ZpzDom<Std64>::init ( Rep& r, const unsigned long long a ) const
{ return r = (Rep)( a >= (uint64_t)_p ? a % (uint64_t)_p : a);
}

inline  ZpzDom<Std64>::Rep&  ZpzDom<Std64>::init ( Rep& r, const double a ) const
{ return init(r, (int64_t)a);
}

inline  ZpzDom<Std64>::Rep&  ZpzDom<Std64>::init ( Rep& r, const float a ) const
{ return init(r, (double)a);
}



inline  ZpzDom<Std64>::Rep&  ZpzDom<Std64>::init ( Rep& r, const long long a ) const
{
  int sign; uint64_t ua;
  if (a <0) { sign =-1; ua = (unsigned int)-a;}
  else { ua = (unsigned int)a; sign =1; }
  r = (Rep) ( (ua >=_p) ? ua % (uint64_t)_p : ua) ;
  if (r && (sign ==-1)) r = (Rep)_p - r;
  return r;
}

/*
inline  ZpzDom<Std64>::Rep&  ZpzDom<Std64>::assign ( Rep& r, const long a ) const
{  return r = (Rep)a;
}

inline  ZpzDom<Std64>::Rep&  ZpzDom<Std64>::assign ( Rep& r, const int a ) const
{ return assign( r, (long)a); }

inline  ZpzDom<Std64>::Rep&  ZpzDom<Std64>::assign ( Rep& r, const unsigned long a ) const
{ return r = (Rep)a; }

inline  ZpzDom<Std64>::Rep&  ZpzDom<Std64>::assign
  ( Rep& r, const unsigned int a ) const
{ return assign(r, (unsigned long)a); }
*/

inline  ZpzDom<Std64>::Rep&  ZpzDom<Std64>::assign
  ( Rep& r, const Rep a ) const
{ return r = a; }

template< class RandIter >
inline  ZpzDom<Std64>::Rep& ZpzDom<Std64>::random(RandIter& g, Rep& a) const {
	return init(a, g());
}

template< class RandIter >
inline  ZpzDom<Std64>::Rep& ZpzDom<Std64>::random(RandIter& g, Rep& a, const Rep& b) const {
	return init(a, g());
}

template< class RandIter >
inline  ZpzDom<Std64>::Rep& ZpzDom<Std64>::random(RandIter& g, Rep& a, long b) const {
	return init(a, g() %(uint64_t) b);
}

template< class RandIter >
inline  ZpzDom<Std64>::Rep& ZpzDom<Std64>::nonzerorandom(RandIter& g, Rep& a) const {
	while (isZero(init(a, g()))) {};
	return a;
}

template< class RandIter >
inline  ZpzDom<Std64>::Rep& ZpzDom<Std64>::nonzerorandom(RandIter& g, Rep& a, const Rep& b) const {
	while (isZero(init(a, g()))) {};
	return a;
}

template< class RandIter >
inline  ZpzDom<Std64>::Rep& ZpzDom<Std64>::nonzerorandom(RandIter& g, Rep& a, long b) const {
	while (isZero(init(a, g() %(uint64_t) b))) {};
	return a;
}

inline void ZpzDom<Std64>::init ( const size_t sz, Array r, constArray a ) const
{
	for ( size_t i = sz; --i; ) {
		r[i] = a[i];
	}
}

inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::init ( Rep& r ) const
{ return r = zero; }

inline void ZpzDom<Std64>::dotprod
  ( Rep& r, const int bound, const size_t sz, constArray a, constArray b ) const
{
  unsigned int stride = 1;
  if ((int64_t)bound < Signed_Trait<Rep>::max() )
   stride = (unsigned int) ( GIVARO_MAXULONG/((unsigned long)bound * (unsigned long)bound) );
  unsigned long dot = (unsigned long) zero; // this is intented !
  if ((sz <10) && (sz <stride)) {
    for(  size_t i= sz; i--; )
#ifdef __x86_64__
      dot += (unsigned long)a[i] * (unsigned long)b[i];
#else
      dot = (unsigned long) (dot + a[i] * b[i]);
#endif
    if (dot > _p) r = (Rep)(dot % (uint64_t)_p);
    else r = (Rep)dot;
    return;
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
    r = (Rep)dot;
    return;
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
  r = (Rep)dot;
}

inline void ZpzDom<Std64>::dotprod
  ( Rep& r, const size_t sz, constArray a, constArray b ) const
{
  ZpzDom<Std64>::dotprod(r, int(_p), sz, a, b);
}


  //  a -> r: int64_t to double
inline void
  ZpzDom<Std64>::i2d ( const size_t sz, double* r, constArray a ) const
{
  for (size_t i=0; i<sz; ++i)  {
	  r[i] = (double) a[i];
  }
}

  //  a -> r: double to int64_t
inline void
  ZpzDom<Std64>::d2i ( const size_t sz, Array r, const double* a ) const
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
inline std::istream& ZpzDom<Std64>::read (std::istream& s)
{
  char ch;
  s >> std::ws >> ch;
  if (ch != '(')
//    GivError::throw_error( GivBadFormat("ZpzDom<Std64>::read: syntax error: no '('"));
    std::cerr << "GivBadFormat(ZpzDom<Std64>::read: syntax error: no '('))" << std::endl;

  s >> std::ws >> ch;
  if (ch != 'z')
//    GivError::throw_error( GivBadFormat("ZpzDom<Std64>::read: bad domain object"));
    std::cerr << "GivBadFormat(ZpzDom<Std64>::read: bad domain object))" << std::endl;

  s >> std::ws >> ch;
  if (ch != ',')
//    GivError::throw_error( GivBadFormat("ZpzDom<Std64>::read: syntax error: no ','"));
    std::cerr << "GivBadFormat(ZpzDom<Std64>::read: syntax error: no ',')) " << std::endl;

  s >> std::ws >> _p;


  s >> std::ws >> ch;
  if (ch != ')')
//    GivError::throw_error( GivBadFormat("ZpzDom<Std64>::read: syntax error: no ')'"));
    std::cerr << "GivBadFormat(ZpzDom<Std64>::read: syntax error: no ')')) " << std::endl;

  return s;
}

inline std::ostream& ZpzDom<Std64>::write (std::ostream& s ) const
{
  return s << "Std64 Givaro Z/pZ modulo " << residu();
}

inline std::istream& ZpzDom<Std64>::read (std::istream& s, Rep& a) const
{
  s >> a;
  init(a, a);
  return s;
}

template <class XXX> inline XXX& ZpzDom<Std64>::convert (XXX& s, const Rep a) const
{
  return s = XXX(a);
}

inline std::ostream& ZpzDom<Std64>::write (std::ostream& s, const Rep a) const
{
  return s << a;
}


inline Integer& ZpzDom<Std64>::write (Integer& r, const Rep a) const
{
  return r = Integer(a);
}

#endif // __DONOTUSE_Givaro_SIXTYFOUR__

} // namespace Givaro

#endif // __GIVARO_zpz64std_INL
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
