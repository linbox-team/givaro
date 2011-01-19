// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givzpzInt.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: JG Dumas
// $Id: givzpzInt.inl,v 1.9 2011-01-19 18:29:09 bboyer Exp $
// ==========================================================================
// Description:

// ---------
// -- normalized operations
// ---------

// r = a*b
// #define __GIVARO_ZPZInteger_N_MUL(r,p,a,b) { r = a*b % p; }
#define __GIVARO_ZPZInteger_N_MUL(r,p,a,b) { r = a; r*=b; r %= p; }
// r *= a
//#define __GIVARO_ZPZInteger_N_MULIN(r,p,a) {  r = (r*a % p);  }
#define __GIVARO_ZPZInteger_N_MULIN(r,p,a) {  r *= a; r %= p;  }

// r = a - b
//#define __GIVARO_ZPZInteger_N_SUB(r,p,a,b) { r = (a-b); r= (r < 0 ? r+p : r); }
#define __GIVARO_ZPZInteger_N_SUB(r,p,a,b) { r = (a-b); if (r < 0 ) r+=p; }
// r -= a
// #define __GIVARO_ZPZInteger_N_SUBIN(r,p,a) { r -= a; r= (r < 0 ? r+p : r); }
#define __GIVARO_ZPZInteger_N_SUBIN(r,p,a) { r -= a; if (r < 0 ) r+=p; }

// r = a+b
// #define __GIVARO_ZPZInteger_N_ADD(r,p,a,b) { r = (a+b); r= (r < p ? r : r-p); }
#define __GIVARO_ZPZInteger_N_ADD(r,p,a,b) { r = (a+b); if (r >= p) r-=p; }
// r += a
// #define __GIVARO_ZPZInteger_N_ADDIN(r,p,a) { r += a;  r= (r < p ? r : r-p); }
#define __GIVARO_ZPZInteger_N_ADDIN(r,p,a) { r += a;  if (r >= p) r-=p; }

// r <- a*b+c % p
// #define __GIVARO_ZPZInteger_N_MULADD(r,p,a,b,c) { r = (a*b+c) % p;  }
#define __GIVARO_ZPZInteger_N_MULADD(r,p,a,b,c) { r = a; r*=b; r+=c; r %= p;  }

// #define __GIVARO_ZPZInteger_N_MULADDIN(r,p,a,b) { r = (a*b+r) % p;  }
#define __GIVARO_ZPZInteger_N_MULADDIN(r,p,a,b) { r += (a*b); r %= p;  }

// a*b-c
//#define __GIVARO_ZPZInteger_N_MULSUB(r,p,a,b,c) { r = (a*b+p-c); r= (r<p ? r : r % p);  }
#define __GIVARO_ZPZInteger_N_MULSUB(r,p,a,b,c) { r = a; r*=b; r+=p; r-=c; if (r>=p) r %= p;  }
// a*b-c
//#define __GIVARO_ZPZInteger_N_SUBMULIN(r,p,a,b) { r -= (a*b); if (r<0) { r+=p; r = (r<0 ? r % p : r); } }
#define __GIVARO_ZPZInteger_N_SUBMULIN(r,p,a,b) { r -= (a*b); if (r<0) { r+=p; if (r<0 ) r %= p ; if (r<0 ) r += p ; } }

#define __GIVARO_ZPZInteger_N_NEG(r,p,a) { r = ( isZero(a) ? zero : p-a); }
#define __GIVARO_ZPZInteger_N_NEGIN(r,p) { r = ( isZero(r) ? zero : p-r); }


inline ZpzDom<Integer>::Residu_t ZpzDom<Integer>::residu( ) const
{ return _p; }



 // ------------------------- Miscellaneous functions

inline int ZpzDom<Integer>::isZero(const Rep& a) const
{ return ::isZero(a); }

inline int ZpzDom<Integer>::isOne(const Rep& a) const
{ return ::isOne(a); }



inline size_t ZpzDom<Integer>::length(const Rep& a) const
{ return ::length(a);}



 // ------------------------- Arithmetic functions




inline ZpzDom<Integer>::Rep& ZpzDom<Integer>::mul (Rep& r, const Rep& a, const Rep& b) const
{
    __GIVARO_ZPZInteger_N_MUL(r,_p,a,b); return r;
}

inline ZpzDom<Integer>::Rep& ZpzDom<Integer>::sub (Rep& r, const Rep& a, const Rep& b) const
{
  __GIVARO_ZPZInteger_N_SUB(r,_p,a,b); return r;
}

inline ZpzDom<Integer>::Rep& ZpzDom<Integer>::add (Rep& r, const Rep& a, const Rep& b) const
{
    __GIVARO_ZPZInteger_N_ADD(r,_p,a,b); return r;
}

inline ZpzDom<Integer>::Rep& ZpzDom<Integer>::neg (Rep& r, const Rep& a) const
{
    __GIVARO_ZPZInteger_N_NEG(r,_p,a); return r;

}

inline ZpzDom<Integer>::Rep& ZpzDom<Integer>::negin (Rep& r) const
{
  __GIVARO_ZPZInteger_N_NEGIN(r,_p);
  return r;
}

inline ZpzDom<Integer>::Rep& ZpzDom<Integer>::inv (Rep& r, const Rep& a) const
{
//  Rep d, v;
//  d = gcd(a, _p, r, v);
//  if (d == -1) negin(r);
//  return r = (r<0)?r + _p:r;
// JGD 03.06.2003
	return ::inv(r,a,_p);
}

inline ZpzDom<Integer>::Rep& ZpzDom<Integer>::div (Rep& r, const Rep& a, const Rep& b) const
{
  Rep ib;
  inv(ib, b);
  __GIVARO_ZPZInteger_N_MUL(r,_p,a,ib);
  return r;
}

 // -- inline array operations between ZpzDom<Integer>::Rep
inline void ZpzDom<Integer>::mul (const size_t sz, Array r, constArray a, constArray b) const
{
  for ( size_t i=sz ; --i ; ) {
    __GIVARO_ZPZInteger_N_MUL(r[i], _p,a[i], b[i]);
  }
}

inline void ZpzDom<Integer>::mul (const size_t sz, Array r, constArray a, const Rep& b) const
{
  for ( size_t i=sz ; --i ; ) {
    __GIVARO_ZPZInteger_N_MUL(r[i], _p, a[i], b);
  }
}

inline void ZpzDom<Integer>::div (const size_t sz, Array r, constArray a, constArray b) const
{
  for ( size_t i=sz ; --i ; ) {
    div( r[i], a[i], b[i]);
  }
}

inline void ZpzDom<Integer>::div (const size_t sz, Array r, constArray a, const Rep& b) const
{
  ZpzDom<Integer>::Rep ib;
  inv(ib, b);
  mul(sz, r, a, ib);
}

inline void ZpzDom<Integer>::add (const size_t sz, Array r, constArray a, constArray b) const
{
  for ( size_t i=sz ; --i ; ) {
    __GIVARO_ZPZInteger_N_ADD(r[i], _p, a[i], b[i]);
  }
}

inline void ZpzDom<Integer>::add (const size_t sz, Array r, constArray a, const Rep& b) const
{
  for ( size_t i=sz ; --i ; ) {
    __GIVARO_ZPZInteger_N_ADD(r[i], _p, a[i], b);
  }
}

inline void ZpzDom<Integer>::sub (const size_t sz, Array r, constArray a, constArray b) const
{
  for ( size_t i=sz ; --i ; ) {
    __GIVARO_ZPZInteger_N_SUB(r[i], _p, a[i], b[i]);
  }
}

inline void ZpzDom<Integer>::sub (const size_t sz, Array r, constArray a, const Rep& b) const
{
  for ( size_t i=sz ; --i ; ) {
    __GIVARO_ZPZInteger_N_SUB(r[i], _p, a[i], b);
  }
}

inline void ZpzDom<Integer>::neg (const size_t sz, Array r, constArray a) const
{
  for ( size_t i=sz ; --i ; ) {
    __GIVARO_ZPZInteger_N_NEG(r[i], _p, a[i]);
  }
}


inline ZpzDom<Integer>::Rep& ZpzDom<Integer>::mulin (Rep& r, const Rep& a) const
{
  __GIVARO_ZPZInteger_N_MULIN(r,_p, a);
  return r;
}

inline ZpzDom<Integer>::Rep& ZpzDom<Integer>::divin (Rep& r, const Rep& a) const
{
  ZpzDom<Integer>::Rep ia;
  inv(ia, a);
  return mulin(r, ia);
}

inline ZpzDom<Integer>::Rep& ZpzDom<Integer>::addin (Rep& r, const Rep& a) const
{
  __GIVARO_ZPZInteger_N_ADDIN(r,_p, a);
  return r;
}

inline ZpzDom<Integer>::Rep& ZpzDom<Integer>::subin (Rep& r, const Rep& a) const
{
  __GIVARO_ZPZInteger_N_SUBIN(r,_p, a);
  return r;
}


inline ZpzDom<Integer>::Rep& ZpzDom<Integer>::invin (Rep& r) const
{
//  Rep d, u, v;
//  d = gcd(r, _p, u, v);
//  if (d == -1) negin(u);
//  return r = (u<0)?u + _p:u;
// JGD 03.06.2003
   Rep t = r;
   return ::inv(r,t,_p);
}

inline ZpzDom<Integer>::Rep& ZpzDom<Integer>::axpy (Rep& r,
						    const Rep& a, const Rep& b, const Rep& c) const
{
  __GIVARO_ZPZInteger_N_MULADD(r, _p, a, b, c);
  return r;
}

inline ZpzDom<Integer>::Rep&  ZpzDom<Integer>::axpyin (Rep& r,
						       const Rep& a, const Rep& b) const
{
  Rep tmp = r;
  __GIVARO_ZPZInteger_N_MULADDIN(tmp, _p, a, b);
  return r = (ZpzDom<Integer>::Rep)tmp;
}


inline void ZpzDom<Integer>::axpy (const size_t sz, Array r,
				   constArray a, constArray x, constArray y) const
{
  for ( size_t i=sz ; --i ; ) {
    __GIVARO_ZPZInteger_N_MULADD(r[i], _p, a[i], x[i], y[i]);
  }
}

inline void ZpzDom<Integer>::axpyin (const size_t sz, Array r,
				     constArray a, constArray x) const
{
  for ( size_t i=sz ; --i ; ) {
    Rep tmp = r[i];
    __GIVARO_ZPZInteger_N_MULADDIN(tmp, _p, a[i], x[i]);
    r[i] = (ZpzDom<Integer>::Rep)tmp;
  }
}

inline ZpzDom<Integer>::Rep&  ZpzDom<Integer>::axmy (Rep& r,
						     const Rep& a, const Rep& b, const Rep& c) const
{
  __GIVARO_ZPZInteger_N_MULSUB(r, _p, a, b, c);
  return r;
}

// r = c - a*b
inline ZpzDom<Integer>::Rep&  ZpzDom<Integer>::maxpy (Rep& r,
						      const Rep& a, const Rep& b, const Rep& c) const
{
  Rep tmp = c;
  __GIVARO_ZPZInteger_N_SUBMULIN(tmp, _p, a, b );
  return r = (ZpzDom<Integer>::Rep)tmp;
}
// r -= a*b
inline ZpzDom<Integer>::Rep&  ZpzDom<Integer>::maxpyin (Rep& r,
						       	const Rep& a, const Rep& b) const
{
  __GIVARO_ZPZInteger_N_SUBMULIN(r, _p, a, b );
  return r;
//   Rep tmp = r;
//   __GIVARO_ZPZInteger_N_SUBMULIN(tmp, _p, a, b );
//   return r = (ZpzDom<Integer>::Rep)tmp;
}
// r = a*b - r
inline ZpzDom<Integer>::Rep&  ZpzDom<Integer>::axmyin (Rep& r,
						       const Rep& a, const Rep& b) const
{
    maxpyin(r,a,b);
    return negin(r);
}


inline void ZpzDom<Integer>::axmy (const size_t sz, Array r,
				   constArray a, constArray x, constArray y) const
{
  for ( size_t i=sz ; --i ; ) {
    __GIVARO_ZPZInteger_N_MULSUB(r[i], _p, a[i], x[i], y[i]);
  }
}

// r -= a*b
inline void ZpzDom<Integer>::maxpyin (const size_t sz, Array r,
				     constArray a, constArray x) const
{
  for ( size_t i=sz ; --i ; ) {
    Rep tmp = r[i];
    __GIVARO_ZPZInteger_N_SUBMULIN(tmp, _p, a[i], x[i]);
    r[i] = (ZpzDom<Integer>::Rep)tmp;
  }
}


// ---------
// -- misc operations
// ---------


inline  ZpzDom<Integer>::Rep&  ZpzDom<Integer>::init ( Rep& r, const double a ) const
{
  int sign; double ua;
  if (a < 0.0) { sign =-1; ua = -a;}
  else { ua = a; sign =1; }
  r = Integer(ua);
  if (r >=_p) r %= _p;
  if (!isZero(r) && (sign == -1)) r = _p - r;
  return r;
}

inline  ZpzDom<Integer>::Rep&  ZpzDom<Integer>::init ( Rep& r, const float a ) const {
    return init(r, (double)a);
}



inline  ZpzDom<Integer>::Rep&  ZpzDom<Integer>::init ( Rep& r, const unsigned long a ) const
{
    r = Integer(a);
    if ( r >= _p ) r %= _p;
    return r ;
}

inline  ZpzDom<Integer>::Rep&  ZpzDom<Integer>::init ( Rep& r, const long a ) const
{
  int sign;
  if (a <0) { sign =-1; r = Integer(-a);}
  else { r = Integer(a); sign =1; }
  if (r >=_p) r %= _p;
  if (!isZero(r) && (sign ==-1)) r = _p - r;
  return r;
}

inline  ZpzDom<Integer>::Rep&  ZpzDom<Integer>::init ( Rep& r, const Integer& a ) const
{
  int sign;
  if (a <0) { sign =-1; r = Integer(-a);}
  else { r = Integer(a); sign =1; }
  if (r >=_p) r %= _p;
  if (!isZero(r) && (sign ==-1)) r = _p - r;
  return r;
}

inline ZpzDom<Integer>::Rep& ZpzDom<Integer>::init( Rep& a, const int i) const { return init(a,(long)i); }

inline ZpzDom<Integer>::Rep& ZpzDom<Integer>::init( Rep& a, const unsigned int i) const { return init(a,(unsigned long)i); }


inline void ZpzDom<Integer>::assign
  ( const size_t sz, Array r, constArray a ) const
{
  for ( size_t i=sz ; --i ; ) {
    if (a[i] <ZpzDom<Integer>::zero) {
       r[i] = a[i] + _p;
       if (r[i] <ZpzDom<Integer>::zero) r[i] %= _p;
    }
    else if (a[i] >_p) {
       r[i] = a[i] - _p;
       if (r[i] >=_p) r[i] %= _p;
    }
    else r[i] = a[i];
  }
}

inline  ZpzDom<Integer>::Rep&  ZpzDom<Integer>::assign ( Rep& r, const long a ) const
{
  return r = Rep(a);
}

inline  ZpzDom<Integer>::Rep&  ZpzDom<Integer>::assign ( Rep& r, const short a ) const
{ return ZpzDom<Integer>::assign( r, (long)a); }

inline  ZpzDom<Integer>::Rep&  ZpzDom<Integer>::assign ( Rep& r, const unsigned long a ) const
{ return r = Rep(a); }

inline  ZpzDom<Integer>::Rep&  ZpzDom<Integer>::assign
  ( Rep& r, const unsigned short a ) const
{ return r = Rep(a); }

inline  ZpzDom<Integer>::Rep&  ZpzDom<Integer>::assign
  ( Rep& r, const Rep& a ) const
{ return r=a; }


inline void ZpzDom<Integer>::init
  ( const size_t sz, Array r, constArray a ) const
{
  for ( size_t i=sz ; --i ; )
       r[i] = a[i];
}

inline ZpzDom<Integer>::Rep& ZpzDom<Integer>::init ( Rep& r ) const
{ return r = zero; }


template< class RandIter >
inline  ZpzDom<Integer>::Rep& ZpzDom<Integer>::random(RandIter& g, Rep& a) const {
	        return init(a, g());
}

template< class RandIter >
inline  ZpzDom<Integer>::Rep& ZpzDom<Integer>::random(RandIter& g, Rep& a, const Rep& b) const {
	        Integer::random(a,b);
                return a %= _p;
}
template< class RandIter >
inline  ZpzDom<Integer>::Rep& ZpzDom<Integer>::random(RandIter& g, Rep& a, long b) const {
	        Integer::random(a,b);
	        return a %= _p;

}

template< class RandIter >
inline  ZpzDom<Integer>::Rep& ZpzDom<Integer>::nonzerorandom(RandIter& g, Rep& a) const {
	        while (isZero( random(g,a) )) {};
		return a;
}

template< class RandIter >
inline  ZpzDom<Integer>::Rep& ZpzDom<Integer>::nonzerorandom(RandIter& g, Rep& a, const Rep& b) const {
	        while (isZero( random(g,a,b))) {};
		return a;
}

template< class RandIter >
inline  ZpzDom<Integer>::Rep& ZpzDom<Integer>::nonzerorandom(RandIter& g, Rep& a, long b) const {
	        while (isZero( random(g,a,b))) {};
		return a;
}


  //  a -> r: int32 to double
inline void
  ZpzDom<Integer>::i2d ( const size_t sz, double* r, constArray a ) const
{
  for (size_t i=0; i<sz; ++i) r[i] = a[i];
}

  //  a -> r: double to int32
inline void
  ZpzDom<Integer>::d2i ( const size_t sz, Array r, const double* a ) const
{
  union d_2_l {
    double d;
    int32 r[2];
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
inline std::istream& ZpzDom<Integer>::read (std::istream& s)
{
  char ch;
  s >> std::ws >> ch;
  if (ch != '(')
//    GivError::throw_error( GivBadFormat("ZpzDom<Integer>::read: syntax error: no '('"));
    std::cerr << "GivBadFormat(ZpzDom<Integer>::read: syntax error: no '('))" << std::endl;

  s >> std::ws >> ch;
  if (ch != 'z')
//    GivError::throw_error( GivBadFormat("ZpzDom<Integer>::read: bad domain object"));
    std::cerr << "GivBadFormat(ZpzDom<Integer>::read: bad domain object))" << std::endl;

  s >> std::ws >> ch;
  if (ch != ',')
//    GivError::throw_error( GivBadFormat("ZpzDom<Integer>::read: syntax error: no ','"));
    std::cerr << "GivBadFormat(ZpzDom<Integer>::read: syntax error: no ',')) " << std::endl;

  s >> std::ws >> _p;

  s >> std::ws >> ch;
  if (ch != ')')
//    GivError::throw_error( GivBadFormat("ZpzDom<Integer>::read: syntax error: no ')'"));
    std::cerr << "GivBadFormat(ZpzDom<Integer>::read: syntax error: no ')')) " << std::endl;

  return s;
}

inline std::ostream& ZpzDom<Integer>::write (std::ostream& s ) const
{
  return s << "Arbitrary length (z," << residu() << ')';
}

inline std::istream& ZpzDom<Integer>::read (std::istream& s, Rep& a) const
{
  s >> a;
  init(a, a);
  return s;
}

inline std::ostream& ZpzDom<Integer>::write (std::ostream& s, const Rep& a) const
{
  return s << a;
}
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
