// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givzpz16table1.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: J.G. Dumas$
// Modified by Pascal Giorgi 2002/04/24
// $Id: givzpz16table1.inl,v 1.12 2011-02-02 17:16:43 bboyer Exp $
// ==========================================================================
// Description:

// ---------
// -- normalized operations
// ---------

#ifndef __GIVARO_zpz16_table1_INL
#define __GIVARO_zpz16_table1_INL

#define __GIVARO_ZPZ16_LOG_MUL(r,p,a,b) \
  {(r)= _tab_mul[(a) + (b)]; }
#define __GIVARO_ZPZ16_LOG_DIV(r,p,a,b) \
  {(r)= _tab_div[(a) - (b)]; }
#define __GIVARO_ZPZ16_LOG_INV(r,p,b)   \
  {(r)= _tab_div[ - (b)]; }
#define __GIVARO_ZPZ16_LOG_SUB(r,p,a,b) \
  {(r)= _tab_mul[(a) + _tab_subone[(b) - (a)] ]; }
#define __GIVARO_ZPZ16_LOG_ADD(r,p,a,b) \
  {(r)=  _tab_mul[(a) + _tab_addone[(b) - (a)] ];}
#define __GIVARO_ZPZ16_LOG_NEG(r,p,a) \
  { r = _tab_neg[(a)];}

/* Pascal Giorgi
   Changing the order of parameters.
*/
#define __GIVARO_ZPZ16_LOG_MULADD(r,p,a,b,c) \
  { __GIVARO_ZPZ16_LOG_MUL(r, p, a, b); __GIVARO_ZPZ16_LOG_ADD(r, p, r, c); }
#define __GIVARO_ZPZ16_LOG_MULSUB(r,p,a,b,c) \
  { __GIVARO_ZPZ16_LOG_MUL(r, p, a, b); __GIVARO_ZPZ16_LOG_SUB(r, p, r, c); }



inline ZpzDom<Log16>::Residu_t ZpzDom<Log16>::residu( ) const
{ return _p; }

inline ZpzDom<Log16>::Rep& ZpzDom<Log16>::mul (Rep& r, const Rep a, const Rep b) const
{
   int32_t tmp;
  __GIVARO_ZPZ16_LOG_MUL(tmp,(int32_t)_p,(int32_t)a,(int32_t)b);
  return r= (ZpzDom<Log16>::Rep)tmp;
}

inline ZpzDom<Log16>::Rep& ZpzDom<Log16>::div (Rep& r, const Rep a, const Rep b) const
{
  __GIVARO_ZPZ16_LOG_DIV(r,_p,a,b);
  return r;
}

inline ZpzDom<Log16>::Rep& ZpzDom<Log16>::sub (Rep& r, const Rep a, const Rep b) const
{
  __GIVARO_ZPZ16_LOG_SUB(r,_p,a,b);
  return r;
}

inline ZpzDom<Log16>::Rep& ZpzDom<Log16>::add (Rep& r, const Rep a, const Rep b) const
{
  __GIVARO_ZPZ16_LOG_ADD(r,_p,a,b);
  return r;
}

inline ZpzDom<Log16>::Rep& ZpzDom<Log16>::neg (Rep& r, const Rep a) const
{
  __GIVARO_ZPZ16_LOG_NEG(r,_p,a);
  return r;
}

inline ZpzDom<Log16>::Rep& ZpzDom<Log16>::inv (Rep& r, const Rep a) const
{
  __GIVARO_ZPZ16_LOG_INV(r,_p,a);
  // GivError::throw_error( GivMathDivZero("Zpz::inv"));
  return r;
}


template< class RandIter >
inline  ZpzDom<Log16>::Rep& ZpzDom<Log16>::random(RandIter& g, Rep& a) const {
	return init(a, g());
}

template< class RandIter >
inline  ZpzDom<Log16>::Rep& ZpzDom<Log16>::random(RandIter& g, Rep& a, const Rep& b) const {
	return init(a, g());
}

template< class RandIter >
inline  ZpzDom<Log16>::Rep& ZpzDom<Log16>::random(RandIter& g, Rep& a, long b) const {
	return init(a, g() %(uint16_t) b);
}

template< class RandIter >
inline  ZpzDom<Log16>::Rep& ZpzDom<Log16>::nonzerorandom(RandIter& g, Rep& a) const {
	while (iszero(init(a, g()))) {};
	return a;
}

template< class RandIter >
inline  ZpzDom<Log16>::Rep& ZpzDom<Log16>::nonzerorandom(RandIter& g, Rep& a, const Rep& b) const {
	while (iszero(init(a, g()))) {};
	return a;
}

template< class RandIter >
inline  ZpzDom<Log16>::Rep& ZpzDom<Log16>::nonzerorandom(RandIter& g, Rep& a, long b) const {
	while (iszero(init(a, g() %(uint16_t) b))) {};
	return a;
}



 // -- inline array operations between ZpzDom<Log16>::Rep
inline void ZpzDom<Log16>::mul (const size_t sz, Array r, constArray a, constArray b) const
{
  for ( size_t i=sz ; --i ; ) {
    __GIVARO_ZPZ16_LOG_MUL(r[i], _p,a[i], b[i]);
  }
}

inline void ZpzDom<Log16>::mul (const size_t sz, Array r, constArray a, Rep b) const
{
  for ( size_t i=sz ; --i ; ) {
    __GIVARO_ZPZ16_LOG_MUL(r[i], _p, a[i], b);
  }
}

inline void ZpzDom<Log16>::div (const size_t sz, Array r, constArray a, constArray b) const
{
  for ( size_t i=sz ; --i ; ) {
    __GIVARO_ZPZ16_LOG_DIV( r[i], _p, a[i], b[i]);
  }
}

inline void ZpzDom<Log16>::div (const size_t sz, Array r, constArray a, Rep b) const
{
  for ( size_t i=sz ; --i ; ) {
    __GIVARO_ZPZ16_LOG_DIV( r[i], _p, a[i], b);
  }
}

inline void ZpzDom<Log16>::add (const size_t sz, Array r, constArray a, constArray b) const
{
  for ( size_t i=sz ; --i ; ) {
    __GIVARO_ZPZ16_LOG_ADD(r[i], _p, a[i], b[i]);
  }
}

inline void ZpzDom<Log16>::add (const size_t sz, Array r, constArray a, Rep b) const
{
  for ( size_t i=sz ; --i ; ) {
    __GIVARO_ZPZ16_LOG_ADD(r[i], _p, a[i], b);
  }
}

inline void ZpzDom<Log16>::sub (const size_t sz, Array r, constArray a, constArray b) const
{
  for ( size_t i=sz ; --i ; ) {
    __GIVARO_ZPZ16_LOG_SUB(r[i], _p, a[i], b[i]);
  }
}

inline void ZpzDom<Log16>::sub (const size_t sz, Array r, constArray a, Rep b) const
{
  for ( size_t i=sz ; --i ; ) {
    __GIVARO_ZPZ16_LOG_SUB(r[i], _p, a[i], b);
  }
}

inline void ZpzDom<Log16>::neg (const size_t sz, Array r, constArray a) const
{
  for ( size_t i=sz ; --i ; ) {
    __GIVARO_ZPZ16_LOG_NEG(r[i], _p, a[i]);
  }
}


inline ZpzDom<Log16>::Rep& ZpzDom<Log16>::mulin (Rep& r, const Rep a) const
{
  __GIVARO_ZPZ16_LOG_MUL(r,_p, r,a);
  return r;
}

inline ZpzDom<Log16>::Rep& ZpzDom<Log16>::divin (Rep& r, const Rep a) const
{
  ZpzDom<Log16>::Rep ia;
  inv(ia, a);
  mulin(r, ia);
  return r;
}

inline ZpzDom<Log16>::Rep& ZpzDom<Log16>::addin (Rep& r, const Rep a) const
{
  __GIVARO_ZPZ16_LOG_ADD(r, _p, r,a);
  return r;
}

inline ZpzDom<Log16>::Rep& ZpzDom<Log16>::subin (Rep& r, const Rep a) const
{
  __GIVARO_ZPZ16_LOG_SUB(r,_p, r,a);
  return r;
}

inline ZpzDom<Log16>::Rep& ZpzDom<Log16>::negin (Rep& r) const
{
  __GIVARO_ZPZ16_LOG_NEG(r,_p,r);
  return r;
}

inline ZpzDom<Log16>::Rep&  ZpzDom<Log16>::invin (Rep& r) const
{
  __GIVARO_ZPZ16_LOG_INV(r,_p,r);
  return r;
}

inline ZpzDom<Log16>::Rep& ZpzDom<Log16>::axpy
 (Rep& r, const Rep a, const Rep b, const Rep c) const
{
  __GIVARO_ZPZ16_LOG_MULADD(r, _p, a, b, c);
  return r;
}

inline ZpzDom<Log16>::Rep& ZpzDom<Log16>::axpyin
 (Rep& r, const Rep a, const Rep b) const
{
  //__GIVARO_ZPZ16_LOG_MULADD(r, _p, a, b, r);
  // Pascal giorgi
  // non consistant because don't perform an axpyin
  // due to the name of parameter
  // the operation performed is 2*a*b

  return axpy(r,a,b,r);
}


inline void ZpzDom<Log16>::axpy
  (const size_t sz, Array r, constArray a, constArray x, constArray y) const
{
  for ( size_t i=sz ; --i ; ) {
    __GIVARO_ZPZ16_LOG_MULADD(r[i], _p, a[i], x[i], y[i]);
  }
}

inline void ZpzDom<Log16>::axpyin
  (const size_t sz, Array r, constArray a, constArray x) const
{
  for ( size_t i=sz ; --i ; ) {
    __GIVARO_ZPZ16_LOG_MULADD(r[i], _p, a[i], x[i], r[i]);
  }
}

// r <- a*b-c
inline void ZpzDom<Log16>::axmy
  (Rep& r, const Rep a, const Rep b, const Rep c) const
{
  __GIVARO_ZPZ16_LOG_MULSUB(r,_p,a,b,c);
}

// r <- r-a*b
inline void ZpzDom<Log16>::maxpyin
  (Rep& r, const Rep a, const Rep b) const
{
    Rep t; __GIVARO_ZPZ16_LOG_MUL(t,_p,a,b);
    this->addin(r,this->negin(t));
}

// r <- c-a*b
inline void ZpzDom<Log16>::maxpy
  (Rep& r, const Rep a, const Rep b, const Rep c) const
{
    Rep t; __GIVARO_ZPZ16_LOG_MUL(t,_p,a,b);
    this->sub(r,c,t);
}


// r <- a*b-r
inline void ZpzDom<Log16>::axmyin (Rep& r,
	   	const Rep a, const Rep b) const
{
    maxpyin(r,a,b);
    negin(r);
}

inline void ZpzDom<Log16>::axmy
  (const size_t sz, Array r, constArray a, constArray x, constArray y) const
{
  for ( size_t i=sz ; --i ; ) {
    __GIVARO_ZPZ16_LOG_MULSUB(r[i], _p, a[i], x[i], y[i]);
  }
}

inline void ZpzDom<Log16>::maxpyin (const size_t sz, Array r,
	   	constArray a, constArray x) const
{
  for ( size_t i=sz ; --i ; ) {
    __GIVARO_ZPZ16_LOG_MULSUB(r[i], _p, a[i], x[i], r[i]);
    __GIVARO_ZPZ16_LOG_NEG(r[i], _p, r[i]);
  }
}

 // ------------------------- Miscellaneous functions

inline int ZpzDom<Log16>::iszero(const Rep a) const
{ return a >= _p; }

inline int ZpzDom<Log16>::isone(const Rep a) const
{ return a == ZpzDom<Log16>::one; }

inline size_t ZpzDom<Log16>::length(const Rep ) const
{ return ZpzDom<Log16>::size_rep;}

inline int ZpzDom<Log16>::isZero( const Rep a ) const {return iszero(a);}
inline int ZpzDom<Log16>::isOne ( const Rep a ) const {return isone(a);}


// ---------
// -- misc operations
// ---------
/*
inline void ZpzDom<Log16>::assign
  ( const size_t sz, Array r, constArray a ) const
{
  for ( size_t i=sz ; --i ; ) {
    if (a[i] <ZpzDom<Log16>::zero) {
       r[i] = a[i] + _p;
       if (r[i] <ZpzDom<Log16>::zero) r[i] = r[i] % _p;
    }
    else if (a[i] >_p) {
       r[i] = a[i] - _p;
       if (r[i] >_p) r[i] = r[i] % _p;
    }
    else r[i] = a[i];
  }
}
*/

inline void ZpzDom<Log16>::assign ( const size_t sz, Array r, constArray a ) const
{
  for ( size_t i=sz ; --i ; )
    r[i] = a[i];
}



// initialized by a degree of the generator.
inline ZpzDom<Log16>::Rep& ZpzDom<Log16>::init ( Rep& r ) const
{ return r = zero; }


inline ZpzDom<Log16>::Rep& ZpzDom<Log16>::assign ( Rep& r, const Rep a ) const
{ return r = a; }

inline ZpzDom<Log16>::Rep& ZpzDom<Log16>::init ( Rep& r, const long a ) const
{
  int sign; unsigned long ua;
  if (a <0) { sign =-1; ua = -a;}
  else { ua = a; sign =1; }
  r = (ua >_p) ? ua % _p : ua;
  if (sign ==-1) r = _p - r;
  return r = _tab_value2rep[r];
}

inline ZpzDom<Log16>::Rep& ZpzDom<Log16>::init ( Rep& r, const int a ) const
{ return ZpzDom<Log16>::init( r, (long)a); }

inline ZpzDom<Log16>::Rep& ZpzDom<Log16>::init ( Rep& r, const unsigned long a ) const
{ r = (a >_p) ? a % _p : a; return r= _tab_value2rep[r];}

inline ZpzDom<Log16>::Rep& ZpzDom<Log16>::init ( Rep& r, const unsigned int a ) const
{ r = (a >_p) ? a % _p : a; return r= _tab_value2rep[r];}

inline ZpzDom<Log16>::Rep& ZpzDom<Log16>::init ( Rep& r, const uint16_t a ) const
{ r = (a >_p) ? a % _p : a; return r= _tab_value2rep[r];}

inline ZpzDom<Log16>::Rep& ZpzDom<Log16>::init ( Rep& r, const int16_t a ) const
{ return ZpzDom<Log16>::init( r, (long)a); }

inline ZpzDom<Log16>::Rep& ZpzDom<Log16>::init( Rep& a, const double i) const {
	  return init(a,(long)i);
}
inline ZpzDom<Log16>::Rep& ZpzDom<Log16>::init( Rep& a, const float i) const {
	  return init(a,(double)i);
}

inline ZpzDom<Log16>::Rep& ZpzDom<Log16>::init ( Rep& r, const Integer& residu ) const
{
  int16_t tr;
  if (residu <0) {
      // -a = b [p]
      // a = p-b [p]
    if ( residu <= (Integer)(-_p) ) tr = int16_t( (-residu) % _p) ;
    else tr = int16_t(-residu);
    if (tr)
      return r = _tab_value2rep[ _p - (uint16_t)tr ];
    else
      return r = zero;
  } else {
    if (residu >= (Integer)_p ) tr =   int16_t(residu % _p) ;
    else tr = int16_t(residu);
    return r = _tab_value2rep[tr];
  }
}



inline void ZpzDom<Log16>::dotprod
  ( Rep& r, const int bound, const size_t sz, constArray a, constArray b ) const
{
  unsigned int stride = 1;
  if ((unsigned long)bound < GIVARO_MAXUINT16)
   stride = GIVARO_MAXUINT32/((unsigned long)bound * (unsigned long)bound);
  uint32_t dot = zero;
  if ((sz <10) && (sz <stride)) {
    for(  int i= sz-1; i>=0; --i)
      dot += _tab_rep2value[a[i]] * _tab_rep2value[b[i]];
    if (dot > _p) r = _tab_value2rep[(Rep)(dot % _p)];
    else r = _tab_value2rep[dot];
    return;
  }
  unsigned int i_begin=0;
  stride &= ~0x1;
  if (stride ==0) {
    for(  int i= sz-1; i>0; --i) {
      dot += _tab_rep2value[a[i]] * _tab_rep2value[b[i]];
      if (dot>_p) dot %= _p;
    }
    r = _tab_value2rep[dot];
    return;
  }
  do {
    size_t min_sz = ((sz-i_begin) < stride ? (sz-i_begin) : stride);
    if ( (min_sz & 0x1) !=0) {
      min_sz--; i_begin++;
      dot += _tab_rep2value[a++[min_sz]] * _tab_rep2value[b++[min_sz]];
    }
    if (min_sz > 1)
      for(  size_t i= min_sz; i>0; --i, --i, ++a, ++a, ++b, ++b )
      {
        dot += _tab_rep2value[a[0]] * _tab_rep2value[b[0]];
        dot += _tab_rep2value[a[1]] * _tab_rep2value[b[1]];
      }
    if (dot>_p) dot %= _p;
    i_begin += min_sz;
  } while (i_begin <sz);
  r = _tab_value2rep[dot];
}

inline void ZpzDom<Log16>::dotprod
  ( Rep& r, const size_t sz, constArray a, constArray b ) const
{
  return ZpzDom<Log16>::dotprod(r, _p, sz, a, b);
}


  //  a -> r: int16_t to double
inline void
  ZpzDom<Log16>::i2d ( const size_t sz, double* r, constArray a ) const
{
  for (size_t i=0; i<sz; ++i) r[i] = _tab_rep2value[a[i]];
}

  //  a -> r: double to int16_t
inline void
  ZpzDom<Log16>::d2i ( const size_t sz, Array r, const double* a ) const
{
  union d_2_l {
    double d;
    int32_t r[2];
  };
  static const double offset = 4503599627370496.0; // 2^52
  size_t i=sz-1;
//warning todo while
//do
label1:
  {
      d_2_l tmp;
      // - normalization: put fractional part at the end of the representation
      tmp.d = a[i] + offset;
      r[i--] = _tab_value2rep[(tmp.r[1] >_p ? tmp.r[1] : tmp.r[1] % _p)];
  }
  // while (i!=0)
  if (i >0) goto label1;
  //for (size_t i=sz-1; i>=0; --i)
}


 // -- Input: (z, <_p>)
inline std::istream& ZpzDom<Log16>::read (std::istream& s)
{
  char ch;
  s >> std::ws >> ch;
//   if (ch != '(')
//     GivError::throw_error( GivBadFormat("ZpzDom<Log16>::read: syntax error: no '('"));
  if (ch != '(')
      std::cerr << "ZpzDom<Log16>::read: syntax error: no '('" << std::endl;

  s >> std::ws >> ch;
//   if (ch != 'z')
//     GivError::throw_error( GivBadFormat("ZpzDom<Log16>::read: bad domain object"));
  if (ch != 'z')
    std::cerr << "ZpzDom<Log16>::read: bad domain object" << std::endl ;

  s >> std::ws >> ch;
//   if (ch != ',')
//     GivError::throw_error( GivBadFormat("ZpzDom<Log16>::read: syntax error: no ','"));
  if (ch != ',')
      std::cerr << "ZpzDom<Log16>::read: syntax error: no ','" << std::endl;


  s >> std::ws >> _p;

  s >> std::ws >> ch;
//   if (ch != ')')
//     GivError::throw_error( GivBadFormat("ZpzDom<Log16>::read: syntax error: no ')'"));
  if (ch != ')')
      std::cerr << "ZpzDom<Log16>::read: syntax error: no ')'" << std::endl;

  return s;
}

inline std::ostream& ZpzDom<Log16>::write (std::ostream& s ) const
{
  return s << "Log16 Givaro Z/pZ modulo " << residu();
}

inline std::istream& ZpzDom<Log16>::read (std::istream& s, Rep& a) const
{
  int tmp; //dpritcha
  s >> tmp;
  tmp %= _p;
  if (tmp < 0) tmp += _p;
  a = _tab_value2rep[tmp];
  return s;
}

inline std::ostream& ZpzDom<Log16>::write (std::ostream& s, const Rep a) const
{
  if (a >= _p) return s << '0';
  return s << _tab_rep2value[a]; //dpritcha
}

#endif // __GIVARO_zpz16_table1_INL
