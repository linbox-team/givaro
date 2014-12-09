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

inline Modular<int64_t>::Residu_t Modular<int64_t>::residu( ) const
{ return _p; }

inline Modular<int64_t>::Rep& Modular<int64_t>::mul (Rep& r, const Rep a, const Rep b) const
{
   int64_t tmp;
  __GIVARO_ZPZ64_N_MUL(tmp,(int64_t)_p,(int64_t)a,(int64_t)b);
  return r = (Modular<int64_t>::Rep)tmp;
}

inline Modular<int64_t>::Rep& Modular<int64_t>::neg (Rep& r, const Rep a) const
{
   int64_t tmp;
  __GIVARO_ZPZ64_N_NEG(tmp,(int64_t)_p,(int64_t)a);
  return r = (Modular<int64_t>::Rep)tmp;
}

inline Modular<int64_t>::Rep& Modular<int64_t>::inv (Rep& r, const Rep a) const
{
//    int64_t d, u, v;
   int64_t u;
  Modular<int64_t>::invext(u, a, (int64_t)_p);
  return r = (u<0)?(Modular<int64_t>::Rep)(u+(int64_t)_p):(Modular<int64_t>::Rep)u;
 }

inline Modular<int64_t>::Rep& Modular<int64_t>::div (Rep& r, const Rep a, const Rep b) const
{
   int64_t tmp;
   int64_t ib;
  inv(ib, b);
  __GIVARO_ZPZ64_N_MUL(tmp,(int64_t)_p,(int64_t)a,(int64_t)ib);
  return r = (Modular<int64_t>::Rep)tmp;
}

inline Modular<int64_t>::Rep& Modular<int64_t>::sub (Rep& r, const Rep a, const Rep b) const
{
   int64_t tmp;
  __GIVARO_ZPZ64_N_SUB(tmp,(int64_t)_p,(int64_t)a,(int64_t)b);
  return r = (Modular<int64_t>::Rep)tmp;
}

inline Modular<int64_t>::Rep& Modular<int64_t>::add (Rep& r, const Rep a, const Rep b) const
{
   int64_t tmp;
  __GIVARO_ZPZ64_N_ADD(tmp,(int64_t)_p,(int64_t)a,(int64_t)b);
  return r = (Modular<int64_t>::Rep)tmp;
}

 // -- inline array operations between Modular<int64_t>::Rep
inline void Modular<int64_t>::mul (const size_t sz, Array r, constArray a, constArray b) const
{
	for ( size_t i = sz; --i; ) {
		int64_t tmp;
		__GIVARO_ZPZ64_N_MUL(tmp, (int64_t)_p,(int64_t)a[i], (int64_t)b[i]);
		r[i] = (Modular<int64_t>::Rep)tmp;
	}
}

inline void Modular<int64_t>::mul (const size_t sz, Array r, constArray a, Rep b) const
{
	for ( size_t i = sz; --i; ) {
		int64_t tmp;
		__GIVARO_ZPZ64_N_MUL(tmp, (int64_t)_p, (int64_t)a[i], (int64_t)b);
		r[i] = (Modular<int64_t>::Rep)tmp;
	}
}

inline void Modular<int64_t>::div (const size_t sz, Array r, constArray a, constArray b) const
{
	for ( size_t i = sz; --i; ) {
		div( r[i], a[i], b[i]);
	}
}

inline void Modular<int64_t>::div (const size_t sz, Array r, constArray a, Rep b) const
{
  Modular<int64_t>::Rep ib;
  inv(ib, b);
  mul(sz, r, a, ib);
}

inline void Modular<int64_t>::add (const size_t sz, Array r, constArray a, constArray b) const
{
	for ( size_t i = sz; --i; ) {
		int64_t tmp;
		__GIVARO_ZPZ64_N_ADD(tmp, (int64_t)_p, (int64_t)a[i], (int64_t)b[i]);
		r[i] = (Modular<int64_t>::Rep)tmp;
	}
}

inline void Modular<int64_t>::add (const size_t sz, Array r, constArray a, Rep b) const
{
	for ( size_t i = sz; --i; ) {
		int64_t tmp;
		__GIVARO_ZPZ64_N_ADD(tmp, (int64_t)_p, (int64_t)a[i], (int64_t)b);
		r[i] = (Modular<int64_t>::Rep)tmp;
	}
}

inline void Modular<int64_t>::sub (const size_t sz, Array r, constArray a, constArray b) const
{
	for ( size_t i = sz; --i; ) {
		int64_t tmp;
		__GIVARO_ZPZ64_N_SUB(tmp, (int64_t)_p, (int64_t)a[i], (int64_t)b[i]);
		r[i] = (Modular<int64_t>::Rep)tmp;
	}
}

inline void Modular<int64_t>::sub (const size_t sz, Array r, constArray a, Rep b) const
{
	for ( size_t i = sz; --i; ) {
		int64_t tmp;
		__GIVARO_ZPZ64_N_SUB(tmp, (int64_t)_p, (int64_t)a[i], (int64_t)b);
		r[i] = (Modular<int64_t>::Rep)tmp;
	}
}

inline void Modular<int64_t>::neg (const size_t sz, Array r, constArray a) const
{
	for ( size_t i = sz; --i; ) {
		int64_t tmp;
		__GIVARO_ZPZ64_N_NEG(tmp, (int64_t)_p, (int64_t)a[i]);
		r[i] = (Modular<int64_t>::Rep)tmp;
	}
}


inline Modular<int64_t>::Rep& Modular<int64_t>::mulin (Rep& r, const Rep a) const
{
   int64_t tmp = (int64_t)r;
  __GIVARO_ZPZ64_N_MULIN(tmp,(int64_t)_p, (int64_t)a);
  return r = (Modular<int64_t>::Rep)tmp;
}

inline Modular<int64_t>::Rep& Modular<int64_t>::divin (Rep& r, const Rep a) const
{
  Modular<int64_t>::Rep ia;
  inv(ia, a);
  return mulin(r, ia);
}

inline Modular<int64_t>::Rep& Modular<int64_t>::addin (Rep& r, const Rep a) const
{
   int64_t tmp = (int64_t)r;
  __GIVARO_ZPZ64_N_ADDIN(tmp,(int64_t)_p, (int64_t)a);
  return r = (Modular<int64_t>::Rep)tmp;
}

inline Modular<int64_t>::Rep& Modular<int64_t>::subin (Rep& r, const Rep a) const
{
   int64_t tmp = (int64_t)r;
  __GIVARO_ZPZ64_N_SUBIN(tmp,(int64_t)_p, (int64_t)a);
  return r = (Modular<int64_t>::Rep)tmp;
}


inline Modular<int64_t>::Rep& Modular<int64_t>::negin (Rep& r) const
{
  __GIVARO_ZPZ64_N_NEGIN(r,(int64_t)_p);
  return r;
}

inline Modular<int64_t>::Rep& Modular<int64_t>::invin (Rep& r) const
{
   int64_t u;
  Modular<int64_t>::invext(u, r, (int64_t)_p);
  return r = (u<0)?(Modular<int64_t>::Rep)(u+(int64_t)_p):(Modular<int64_t>::Rep)u;
}


inline Modular<int64_t>::Rep& Modular<int64_t>::axpy
 (Rep& r, const Rep a, const Rep b, const Rep c) const
{
   int64_t tmp;
  __GIVARO_ZPZ64_N_MULADD(tmp, (int64_t)_p, (int64_t)a, (int64_t)b, (int64_t)c);
  return r = (Modular<int64_t>::Rep)tmp;
}

inline Modular<int64_t>::Rep& Modular<int64_t>::axpyin
 (Rep& r, const Rep a, const Rep b) const
{
   int64_t tmp = (int64_t)r;
  __GIVARO_ZPZ64_N_MULADDIN(tmp, (int64_t)_p, (int64_t)a, (int64_t)b);
  return r = (Modular<int64_t>::Rep)tmp;
}


inline void Modular<int64_t>::axpy
  (const size_t sz, Array r, constArray a, constArray x, constArray y) const
{
	for ( size_t i = sz; --i; ) {
		int64_t tmp;
		__GIVARO_ZPZ64_N_MULADD(tmp, (int64_t)_p, (int64_t)a[i], (int64_t)x[i], (int64_t)y[i]);
		r[i] = (Modular<int64_t>::Rep)tmp;
	}
}

inline void Modular<int64_t>::axpyin
  (const size_t sz, Array r, constArray a, constArray x) const
{
	for ( size_t i = sz; --i; ) {
		int64_t tmp = (int64_t)r[i];
		__GIVARO_ZPZ64_N_MULADDIN(tmp, (int64_t)_p, (int64_t)a[i], (int64_t)x[i]);
		r[i] = (Modular<int64_t>::Rep)tmp;
	}
}

inline Modular<int64_t>::Rep&  Modular<int64_t>::axmy
 (Rep& r, const Rep a, const Rep b, const Rep c) const
{
   int64_t tmp;
  __GIVARO_ZPZ64_N_MULSUB(tmp, (int64_t)_p, (int64_t)a, (int64_t)b, (int64_t)c);
  return r = (Modular<int64_t>::Rep)tmp;
}

// r = c - a*b
inline Modular<int64_t>::Rep&  Modular<int64_t>::maxpy
 (Rep& r, const Rep a, const Rep b, const Rep c) const
{
   int64_t tmp = (int64_t)c;
  __GIVARO_ZPZ64_N_SUBMULIN(tmp, (int64_t)_p, (int64_t)a, (int64_t)b );
  return r = (Modular<int64_t>::Rep)tmp;
}

#if 0
inline Modular<int64_t>::Rep& Modular<int64_t>::maxpy (Rep& r,
						 const Rep a, const Rep b, const Rep c) const
{
	int64_t tmp;
	__GIVARO_ZPZ64_N_MUL(tmp, (int64_t)_p, (int64_t)a, (int64_t)b);
	__GIVARO_ZPZ64_N_SUB(r, (int64_t)_p, (int64_t)c, tmp);
	return r;
}
#endif

// r -= a*b
inline Modular<int64_t>::Rep&  Modular<int64_t>::maxpyin (Rep& r, const Rep a, const Rep b) const
{
  __GIVARO_ZPZ64_N_SUBMULIN(r, (int64_t)_p, (int64_t)a, (int64_t)b );
  return r;
//    int64_t tmp = (int64_t)r;
//   __GIVARO_ZPZ64_N_SUBMULIN(tmp, (int64_t)_p, (int64_t)a, (int64_t)b );
//   return r = (Modular<int64_t>::Rep)tmp;
}

// r = a*b - r
inline Modular<int64_t>::Rep&  Modular<int64_t>::axmyin (Rep& r, const Rep a, const Rep b) const
{
    maxpyin(r,a,b);
    return negin(r);
}


inline void Modular<int64_t>::axmy (const size_t sz, Array r, constArray a, constArray x, constArray y) const
{
  for ( size_t i = sz; --i; ) {
     int64_t tmp;
    __GIVARO_ZPZ64_N_MULSUB(tmp, (int64_t)_p, (int64_t)a[i], (int64_t)x[i], (int64_t)y[i]);
    r[i] = (Modular<int64_t>::Rep)tmp;
  }
}

// r = a*b - r
inline void Modular<int64_t>::maxpyin (const size_t sz, Array r, constArray a, constArray x) const
{
	for ( size_t i = sz; --i; ) {
		int64_t tmp = (int64_t)r[i];
		__GIVARO_ZPZ64_N_SUBMULIN(tmp, (int64_t)_p, (int64_t)a[i], (int64_t)x[i]);
		r[i] = (Modular<int64_t>::Rep)tmp;
	}
}

 // ------------------------- Miscellaneous functions

inline int Modular<int64_t>::areEqual(const Rep a, const Rep b) const
{ return a == b; }

inline int Modular<int64_t>::areNEqual(const Rep a, const Rep b) const
{ return a != b; }

inline int Modular<int64_t>::isZero(const Rep a) const
{ return a == Modular<int64_t>::zero; }

inline int Modular<int64_t>::isnzero(const Rep a) const
{ return a != Modular<int64_t>::zero; }

inline int Modular<int64_t>::isOne(const Rep a) const
{ return a == Modular<int64_t>::one; }

inline int Modular<int64_t>::isMOne(const Rep a) const
{ return a == Modular<int64_t>::mOne; }


inline size_t Modular<int64_t>::length(const Rep ) const
{ return Modular<int64_t>::size_rep;}


// ---------
// -- misc operations
// ---------

inline void Modular<int64_t>::assign
  ( const size_t sz, Array r, constArray a ) const
{
    for ( size_t i = sz; --i; )
        r[i] = a[i];
}

inline  Modular<int64_t>::Rep&  Modular<int64_t>::init ( Rep& r, const unsigned long a ) const
{ return r = (Rep)( a >= (uint64_t)_p ? a % (uint64_t)_p : a);
}


inline  Modular<int64_t>::Rep&  Modular<int64_t>::init ( Rep& r, const long a ) const
{
  int64_t sign; uint64_t ua;
  if (a <0) { sign =-1; ua = (unsigned int)-a;}
  else { ua = (unsigned int)a; sign =1; }
  r = (Rep)((ua >=_p) ? ua % (uint64_t)_p : ua);
  if (r && (sign ==-1)) r = (Rep)_p - r;
  return r;
}


inline Modular<int64_t>::Rep&  Modular<int64_t>::init ( Rep& r, const Integer& Residu ) const
{
  if (Residu <0) {
      int64_t tr;
          // -a = b [p]
          // a = p-b [p]
      if ( (-Residu) >= (Integer)(_p) ) tr = int64_t( (-Residu) % (Integer)_p) ;
      else tr = int64_t(-Residu);
      if (tr)
          return r = (Rep)( (uint64_t)_p - (uint64_t)tr ) ;
      else
          return r = zero;
  } else {
      Integer ip(_p);
      if (Residu >= ip ) return r =   int64_t(Residu % ip) ;
      else return r = int64_t(Residu);
  }
}




inline  Modular<int64_t>::Rep& Modular<int64_t>::init( Rep& a, const int i) const { return init(a,(long)i); }
inline  Modular<int64_t>::Rep& Modular<int64_t>::init( Rep& a, const unsigned int i) const { return init(a,(unsigned long)i); }


inline  Modular<int64_t>::Rep&  Modular<int64_t>::init ( Rep& r, const unsigned long long a ) const
{ return r = (Rep)( a >= (uint64_t)_p ? a % (uint64_t)_p : a);
}

inline  Modular<int64_t>::Rep&  Modular<int64_t>::init ( Rep& r, const double a ) const
{ return init(r, (int64_t)a);
}

inline  Modular<int64_t>::Rep&  Modular<int64_t>::init ( Rep& r, const float a ) const
{ return init(r, (double)a);
}



inline  Modular<int64_t>::Rep&  Modular<int64_t>::init ( Rep& r, const long long a ) const
{
  int sign; uint64_t ua;
  if (a <0) { sign =-1; ua = (unsigned int)-a;}
  else { ua = (unsigned int)a; sign =1; }
  r = (Rep) ( (ua >=_p) ? ua % (uint64_t)_p : ua) ;
  if (r && (sign ==-1)) r = (Rep)_p - r;
  return r;
}

/*
inline  Modular<int64_t>::Rep&  Modular<int64_t>::assign ( Rep& r, const long a ) const
{  return r = (Rep)a;
}

inline  Modular<int64_t>::Rep&  Modular<int64_t>::assign ( Rep& r, const int a ) const
{ return assign( r, (long)a); }

inline  Modular<int64_t>::Rep&  Modular<int64_t>::assign ( Rep& r, const unsigned long a ) const
{ return r = (Rep)a; }

inline  Modular<int64_t>::Rep&  Modular<int64_t>::assign
  ( Rep& r, const unsigned int a ) const
{ return assign(r, (unsigned long)a); }
*/

inline  Modular<int64_t>::Rep&  Modular<int64_t>::assign
  ( Rep& r, const Rep a ) const
{ return r = a; }

template< class RandIter >
inline  Modular<int64_t>::Rep& Modular<int64_t>::random(RandIter& g, Rep& a) const {
	return init(a, g());
}

template< class RandIter >
inline  Modular<int64_t>::Rep& Modular<int64_t>::random(RandIter& g, Rep& a, const Rep& ) const {
	return init(a, g());
}

template< class RandIter >
inline  Modular<int64_t>::Rep& Modular<int64_t>::random(RandIter& g, Rep& a, long b) const {
	return init(a, g() %(uint64_t) b);
}

template< class RandIter >
inline  Modular<int64_t>::Rep& Modular<int64_t>::nonzerorandom(RandIter& g, Rep& a) const {
	while (isZero(init(a, g()))) {};
	return a;
}

template< class RandIter >
inline  Modular<int64_t>::Rep& Modular<int64_t>::nonzerorandom(RandIter& g, Rep& a, const Rep& ) const {
	while (isZero(init(a, g()))) {};
	return a;
}

template< class RandIter >
inline  Modular<int64_t>::Rep& Modular<int64_t>::nonzerorandom(RandIter& g, Rep& a, long b) const {
	while (isZero(init(a, g() %(uint64_t) b))) {};
	return a;
}

inline void Modular<int64_t>::init ( const size_t sz, Array r, constArray a ) const
{
	for ( size_t i = sz; --i; ) {
		r[i] = a[i];
	}
}

inline Modular<int64_t>::Rep& Modular<int64_t>::init ( Rep& r ) const
{ return r = zero; }

inline void Modular<int64_t>::dotprod
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

inline void Modular<int64_t>::dotprod
  ( Rep& r, const size_t sz, constArray a, constArray b ) const
{
  Modular<int64_t>::dotprod(r, int(_p), sz, a, b);
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
  return s << "int64_t Givaro Z/pZ modulo " << residu();
}

inline std::istream& Modular<int64_t>::read (std::istream& s, Rep& a) const
{
  Integer tmp;
  s >> tmp;
  init(a, tmp);
  return s;
}

template <class XXX> inline XXX& Modular<int64_t>::convert (XXX& s, const Rep a) const
{
  return s = XXX(a);
}

inline std::ostream& Modular<int64_t>::write (std::ostream& s, const Rep a) const
{
  return s << a;
}


inline Integer& Modular<int64_t>::write (Integer& r, const Rep a) const
{
  return r = Integer(a);
}


} // namespace Givaro

#endif // __DONOTUSE_Givaro_SIXTYFOUR__

#endif // __GIVARO_zpz64std_INL
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
