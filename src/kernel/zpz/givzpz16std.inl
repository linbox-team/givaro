// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givzpz16std.inl,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givzpz16std.inl,v 1.3 2004-07-20 12:03:46 giorgi Exp $
// ==========================================================================
// Description:

// ---------
// -- normalized operations
// ---------


// r = a*b
#define __GIVARO_ZPZ16_N_MUL(r,p,a,b) { r = (a*b); r= (r>=p ? r % (uint16)p : r); }
// r *= a
#define __GIVARO_ZPZ16_N_MULIN(r,p,a) { r *= a; r= (r<p ? r : r % (uint16)p);}

// r = a - b
//#define __GIVARO_ZPZ16_N_SUB(r,p,a,b) { r = (a-b); r= (r < 0 ? r+p : r);}
#define __GIVARO_ZPZ16_N_SUB(r,p,a,b) ( r = a>b? a-b: (p-b)+a )

// r -= a
#define __GIVARO_ZPZ16_N_SUBIN(r,p,a) { r -= a; r= (r < 0 ? r+p : r);}

// r = a+b
#define __GIVARO_ZPZ16_N_ADD(r,p,a,b) { r = (a+b); r= (r < p ? r : r-p);}
// r += a
#define __GIVARO_ZPZ16_N_ADDIN(r,p,a) { r += a;  r= (r < p ? r : r-p);}

// r <- a*b+c % p
#define __GIVARO_ZPZ16_N_MULADD(r,p,a,b,c) \
{ r = (a*b+c); r= (r<p ? r : r % (uint16)p); }

#define __GIVARO_ZPZ16_N_MULADDIN(r,p,a,b) \
{ r += (a*b); r= (r<p ? r : r % (uint16)p);}

// // a*b-c
#define __GIVARO_ZPZ16_N_MULSUB(r,p,a,b,c) \
{ r = (a*b-c); r= (r<p ? (r>0? r : r% (uint16)p) : r % (uint16)p); }
// a*b-c
#define __GIVARO_ZPZ16_N_SUBMULIN(r,p,a,b) \
{ r -= (a*b); if (r<0) { r+=p; r = (r<0 ? r % (uint16)p : r); } }




#define __GIVARO_ZPZ16_N_NEG(r,p,a) { r = (a == 0 ? 0 : p-a); }
#define __GIVARO_ZPZ16_N_NEGIN(r,p) { r = (r == 0 ? 0 : p-r); }


inline ZpzDom<Std16>::ZpzDom<Std16>( )
 : zero(0), one(1), _p(0)
{}

inline ZpzDom<Std16>::ZpzDom<Std16>( Residu_t p )
 : zero(0), one(1), _p(p)
{}

inline ZpzDom<Std16>::ZpzDom<Std16>( const ZpzDom<Std16>& F)
 : zero(0), one(1), _p(F._p)
{}

inline ZpzDom<Std16>::Residu_t ZpzDom<Std16>::residu( ) const
{ return _p; }

inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::mul (Rep& r, const Rep a, const Rep b) const
{ 
  register int32 tmp; 
  __GIVARO_ZPZ16_N_MUL(tmp,(int32)_p,(int32)a,(int32)b); 
  return r = (ZpzDom<Std16>::Rep)tmp; 
}

inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::neg (Rep& r, const Rep a) const
{ 
  register int32 tmp; 
  __GIVARO_ZPZ16_N_NEG(tmp,(int32)_p,(int32)a); 
  return r = (ZpzDom<Std16>::Rep)tmp; 
}

inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::inv (Rep& r, const Rep a) const
{ 
  register int32 u; 
  ZpzDom<Std16>::invext(u, a, _p);
//   if ((d != 1) && (d != -1)) std::cerr << "GivMathDivZero(Zpz::inv)" << std::endl;
  return r = (u<0)?(ZpzDom<Std16>::Rep)(u+_p):(ZpzDom<Std16>::Rep)u; 
}

inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::div (Rep& r, const Rep a, const Rep b) const
{ 
  register int32 tmp; 
  register int16 ib; 
  inv(ib, b);
  __GIVARO_ZPZ16_N_MUL(tmp,(int32)_p,(int32)a,(int32)ib); 
  return r = (ZpzDom<Std16>::Rep)tmp; 
}

inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::sub (Rep& r, const Rep a, const Rep b) const
{ 
  return __GIVARO_ZPZ16_N_SUB(r,(int32)_p,(int32)a,(int32)b); 
}

inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::add (Rep& r, const Rep a, const Rep b) const
{ 
  register int32 tmp; 
  __GIVARO_ZPZ16_N_ADD(tmp,(int32)_p,(int32)a,(int32)b); 
  return r = (ZpzDom<Std16>::Rep)tmp; 
}


 // -- inline array operations between ZpzDom<Std16>::Rep
inline void ZpzDom<Std16>::mul (const size_t sz, Array r, constArray a, constArray b) const
{
  for (register size_t i=sz-1; i!=0; --i) {
    register int32 tmp;
    __GIVARO_ZPZ16_N_MUL(tmp, (int32)_p,(int32)a[i], (int32)b[i]); 
    r[i] = (ZpzDom<Std16>::Rep)tmp;
  }
}

inline void ZpzDom<Std16>::mul (const size_t sz, Array r, constArray a, Rep b) const
{
  for (register size_t i=sz-1; i!=0; --i) {
    register int32 tmp;
    __GIVARO_ZPZ16_N_MUL(tmp, (int32)_p, (int32)a[i], (int32)b);
    r[i] = (ZpzDom<Std16>::Rep)tmp;
  }
}

inline void ZpzDom<Std16>::div (const size_t sz, Array r, constArray a, constArray b) const
{
  for (register size_t i=sz-1; i!=0; --i) {
    div( r[i], a[i], b[i]); 
  }
}

inline void ZpzDom<Std16>::div (const size_t sz, Array r, constArray a, Rep b) const
{
  ZpzDom<Std16>::Rep ib;
  inv(ib, b);
  mul(sz, r, a, ib);
}

inline void ZpzDom<Std16>::add (const size_t sz, Array r, constArray a, constArray b) const
{
  for (register size_t i=sz-1; i!=0; --i) {
    register int32 tmp;
    __GIVARO_ZPZ16_N_ADD(tmp, (int32)_p, (int32)a[i], (int32)b[i]);
    r[i] = (ZpzDom<Std16>::Rep)tmp;
  }
}

inline void ZpzDom<Std16>::add (const size_t sz, Array r, constArray a, Rep b) const
{
  for (register size_t i=sz-1; i!=0; --i) {
    register int32 tmp;
    __GIVARO_ZPZ16_N_ADD(tmp, (int32)_p, (int32)a[i], (int32)b);
    r[i] = (ZpzDom<Std16>::Rep)tmp;
  }
}

inline void ZpzDom<Std16>::sub (const size_t sz, Array r, constArray a, constArray b) const
{
  for (register size_t i=sz-1; i!=0; --i) {
    register int32 tmp;
    __GIVARO_ZPZ16_N_SUB(tmp, (int32)_p, (int32)a[i], (int32)b[i]);
    r[i] = (ZpzDom<Std16>::Rep)tmp;
  }
}

inline void ZpzDom<Std16>::sub (const size_t sz, Array r, constArray a, Rep b) const
{
  for (register size_t i=sz-1; i!=0; --i) {
    register int32 tmp;
    __GIVARO_ZPZ16_N_SUB(tmp, (int32)_p, (int32)a[i], (int32)b);
    r[i] = (ZpzDom<Std16>::Rep)tmp;
  }
}

inline void ZpzDom<Std16>::neg (const size_t sz, Array r, constArray a) const
{
  for (register size_t i=sz-1; i!=0; --i) {
    register int32 tmp;
    __GIVARO_ZPZ16_N_NEG(tmp, (int32)_p, (int32)a[i]);
    r[i] = (ZpzDom<Std16>::Rep)tmp;
  }
}


inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::mulin (Rep& r, const Rep a) const
{ 
  register int32 tmp = (int32)r; 
  __GIVARO_ZPZ16_N_MULIN(tmp,(int32)_p, (int32)a); 
  return r = (ZpzDom<Std16>::Rep)tmp; 
}

inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::divin (Rep& r, const Rep a) const
{ 
  ZpzDom<Std16>::Rep ia;
  inv(ia, a);
  return mulin(r, ia);
}

inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::addin (Rep& r, const Rep a) const
{ 
  register int32 tmp = (int32)r; 
  __GIVARO_ZPZ16_N_ADDIN(tmp,(int32)_p, (int32)a); 
  return r = (ZpzDom<Std16>::Rep)tmp; 
}

inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::subin (Rep& r, const Rep a) const
{ 
  register int32 tmp = (int32)r; 
  __GIVARO_ZPZ16_N_SUBIN(tmp,(int32)_p, (int32)a); 
  return r = (ZpzDom<Std16>::Rep)tmp; 
}


inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::negin (Rep& r) const
{ 
  __GIVARO_ZPZ16_N_NEGIN(r,(int32)_p); 
  return r; 
}

inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::invin (Rep& r) const
{ 
  register int32 u; 
  ZpzDom<Std16>::invext(u, r, _p);
//   if ((d != 1) && (d != -1)) std::cerr << "GivMathDivZero(Zpz::inv)" << std::endl;
  return r = (u<0)?(ZpzDom<Std16>::Rep)(u+_p):(ZpzDom<Std16>::Rep)u; 
}


inline void ZpzDom<Std16>::axpy 
 (Rep& r, const Rep a, const Rep b, const Rep c) const
{ 
  register int32 tmp; 
  __GIVARO_ZPZ16_N_MULADD(tmp, (int32)_p, (int32)a, (int32)b, (int32)c); 
  r = (ZpzDom<Std16>::Rep)tmp; 
}

inline void ZpzDom<Std16>::axpyin 
 (Rep& r, const Rep a, const Rep b) const
{ 
  register int32 tmp = (int32)r; 
  __GIVARO_ZPZ16_N_MULADDIN(tmp, (int32)_p, (int32)a, (int32)b); 
  r = (ZpzDom<Std16>::Rep)tmp; 
}


inline void ZpzDom<Std16>::axpy 
  (const size_t sz, Array r, constArray a, constArray x, constArray y) const
{
  for (register size_t i=sz-1; i!=0; --i) {
    register int32 tmp;
    __GIVARO_ZPZ16_N_MULADD(tmp, (int32)_p, (int32)a[i], (int32)x[i], (int32)y[i]);
    r[i] = (ZpzDom<Std16>::Rep)tmp;
  }
}

inline void ZpzDom<Std16>::axpyin 
  (const size_t sz, Array r, constArray a, constArray x) const
{
  for (register size_t i=sz-1; i!=0; --i) {
    register int32 tmp = (int32)r[i];
    __GIVARO_ZPZ16_N_MULADDIN(tmp, (int32)_p, (int32)a[i], (int32)x[i]);
    r[i] = (ZpzDom<Std16>::Rep)tmp;
  }
}

inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::amxy
 (Rep& r, const Rep a, const Rep b, const Rep c) const
{
  register int32 tmp;
  __GIVARO_ZPZ16_N_MUL(tmp, (int32)_p, (int32)a, (int32)b);
  __GIVARO_ZPZ16_N_SUB(r, (int32)_p, (int32)c, tmp);
  return r;
}

inline void ZpzDom<Std16>::axmy
 (Rep& r, const Rep a, const Rep b, const Rep c) const
{
  register int32 tmp;
  __GIVARO_ZPZ16_N_MULSUB(tmp, (int32)_p, (int32)a, (int32)b, (int32)c);
  r = (ZpzDom<Std16>::Rep)tmp;
}

// r -= a*b
inline void ZpzDom<Std16>::axmyin 
 (Rep& r, const Rep a, const Rep b) const
{
  register int32 tmp = (int32)r;
  __GIVARO_ZPZ16_N_SUBMULIN(tmp, (int32)_p, (int32)a, (int32)b );
  r = (ZpzDom<Std16>::Rep)tmp;
}


inline void ZpzDom<Std16>::axmy
  (const size_t sz, Array r, constArray a, constArray x, constArray y) const
{
  for (register size_t i=sz-1; i!=0; --i) {
    register int32 tmp;
    __GIVARO_ZPZ16_N_MULSUB(tmp, (int32)_p, (int32)a[i], (int32)x[i], (int32)y[i]);
    r[i] = (ZpzDom<Std16>::Rep)tmp;
  }
}

// r -= a*b
inline void ZpzDom<Std16>::axmyin
  (const size_t sz, Array r, constArray a, constArray x) const
{
  for (register size_t i=sz-1; i!=0; --i) {
    register int32 tmp = (int32)r[i];
    __GIVARO_ZPZ16_N_SUBMULIN(tmp, (int32)_p, (int32)a[i], (int32)x[i]);
    r[i] = (ZpzDom<Std16>::Rep)tmp;
  }
}

 // ------------------------- Miscellaneous functions

inline int ZpzDom<Std16>::areEqual(const Rep a, const Rep b) const
{ return a == b; }

inline int ZpzDom<Std16>::areNEqual(const Rep a, const Rep b) const
{ return a != b; }

inline int ZpzDom<Std16>::iszero(const Rep a) const
{ return a == ZpzDom<Std16>::zero; }

inline int ZpzDom<Std16>::isnzero(const Rep a) const
{ return a != ZpzDom<Std16>::zero; }

inline int ZpzDom<Std16>::isone(const Rep a) const
{ return a == ZpzDom<Std16>::one; }

inline int ZpzDom<Std16>::isZero( const Rep a ) const { return iszero(a);}
inline int ZpzDom<Std16>::isOne( const Rep a ) const {return isone(a);}


inline size_t ZpzDom<Std16>::length(const Rep a) const
{ return ZpzDom<Std16>::size_rep;}

// ---------
// -- misc operations
// ---------
inline  ZpzDom<Std16>::Rep&  ZpzDom<Std16>::init ( Rep& r, const unsigned long a ) const
{ return r = (Rep)( a >= (unsigned long)_p ? a % (unsigned long)_p : a); 
}

inline  ZpzDom<Std16>::Rep&  ZpzDom<Std16>::init ( Rep& r, const long a ) const
{
  int sign; long ua;
  if (a <0) { sign =-1; ua = -a;}
  else { ua = a; sign =1; }
  r = (ua >=_p) ? ua % (uint16)_p : ua;
  if (r && (sign ==-1)) r = _p - r;
  return r;
}

inline  ZpzDom<Std16>::Rep& ZpzDom<Std16>::init( Rep& a, const int i) const { return init(a,(long)i); }
inline  ZpzDom<Std16>::Rep& ZpzDom<Std16>::init( Rep& a, const double i) const { return init(a,(long)i); }
inline  ZpzDom<Std16>::Rep& ZpzDom<Std16>::init( Rep& a, const float i) const { return init(a,(double)i); }
inline  ZpzDom<Std16>::Rep& ZpzDom<Std16>::init( Rep& a, const unsigned int i) const { return init(a,(unsigned long)i); }


inline void ZpzDom<Std16>::assign 
  ( const size_t sz, Array r, constArray a ) const
{
  for (register size_t i=sz-1; i!=0; --i) {
    if (a[i] <ZpzDom<Std16>::zero) {
       r[i] = a[i] + _p;
       if (r[i] <ZpzDom<Std16>::zero) r[i] = r[i] % (uint16)_p;
    }
    else if (a[i] >_p) {
       r[i] = a[i] - _p;
       if (r[i] >=_p) r[i] = r[i] % (uint16)_p;
    }
    else r[i] = a[i];
  }
}

inline  ZpzDom<Std16>::Rep&  ZpzDom<Std16>::assign ( Rep& r, const long a ) const
{  return init(r, a);
}


inline  ZpzDom<Std16>::Rep&  ZpzDom<Std16>::assign ( Rep& r, const int a ) const
{ return ZpzDom<Std16>::assign( r, (long)a); }

inline  ZpzDom<Std16>::Rep&  ZpzDom<Std16>::assign ( Rep& r, const unsigned long a ) const
{ return init(r,a); }

inline  ZpzDom<Std16>::Rep&  ZpzDom<Std16>::assign 
  ( Rep& r, const unsigned int a ) const
{ return init(r, (unsigned long)a); }

inline  ZpzDom<Std16>::Rep&  ZpzDom<Std16>::assign 
  ( Rep& r, const Rep a ) const
{ return assign(r, (long)a); }

template< class RandIter >
inline  ZpzDom<Std16>::Rep& ZpzDom<Std16>::random(RandIter& g, Rep& a) const {
	return init(a, g());
}

template< class RandIter >
inline  ZpzDom<Std16>::Rep& ZpzDom<Std16>::random(RandIter& g, Rep& a, const Rep& b) const {
	return init(a, g());
}

template< class RandIter >
inline  ZpzDom<Std16>::Rep& ZpzDom<Std16>::random(RandIter& g, Rep& a, long b) const {
	return init(a, g() %(uint16) b);
}

template< class RandIter >
inline  ZpzDom<Std16>::Rep& ZpzDom<Std16>::nonzerorandom(RandIter& g, Rep& a) const {
	while (iszero(init(a, g()))) {};
	return a;
}

template< class RandIter >
inline  ZpzDom<Std16>::Rep& ZpzDom<Std16>::nonzerorandom(RandIter& g, Rep& a, const Rep& b) const {
	while (iszero(init(a, g()))) {};
	return a;
}

template< class RandIter >
inline  ZpzDom<Std16>::Rep& ZpzDom<Std16>::nonzerorandom(RandIter& g, Rep& a, long b) const {
	while (iszero(init(a, g() %(uint16) b))) {};
	return a;
}

inline void ZpzDom<Std16>::init 
  ( const size_t sz, Array r, constArray a ) const
{
  for (register size_t i=sz-1; i!=0; --i) 
       r[i] = a[i];
}

inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::init ( Rep& r ) const
{ return r = zero; }

inline void ZpzDom<Std16>::dotprod 
  ( Rep& r, const int bound, const size_t sz, constArray a, constArray b ) const
{
  unsigned int stride = 1;
  if ((unsigned long)bound < GIVARO_MAXUINT16)
   stride = GIVARO_MAXULONG/((unsigned long)bound * (unsigned long)bound);
  unsigned long dot = zero;
  if ((sz <10) && (sz <stride)) {
    for( register int i= sz-1; i>=0; --i) 
      dot += a[i] * b[i]; 
    if (dot > _p) r = (Rep)(dot % (uint16)_p);
    else r = (Rep)dot;
    return;
  }
  unsigned int i_begin=0;
  stride &= ~0x1;
  if (stride ==0) {
    for( register int i= sz-1; i>0; --i) {
      dot += a[i] * b[i];
      if (dot>_p) dot %= _p;
    }
    r = (Rep)dot;
    return;
  }
  do {
    size_t min_sz = ((sz-i_begin) < stride ? (sz-i_begin) : stride);
    if (min_sz & 0x1 !=0) 
      { min_sz--; i_begin++; dot += a++[min_sz] * b++[min_sz]; }
    if (min_sz > 1) 
      for( register size_t i= min_sz; i>0; --i, --i, ++a, ++a, ++b, ++b ) 
      {
        dot += a[0] * b[0]; 
        dot += a[1] * b[1];
      }
    if (dot>_p) dot %= (uint16)_p;
    i_begin += min_sz;
  } while (i_begin <sz);
  r = (Rep)dot;
}

inline void ZpzDom<Std16>::dotprod
  ( Rep& r, const size_t sz, constArray a, constArray b ) const
{
  return ZpzDom<Std16>::dotprod(r, _p, sz, a, b);
}


  //  a -> r: int16 to double
inline void 
  ZpzDom<Std16>::i2d ( const size_t sz, double* r, constArray a ) const
{
  for (size_t i=0; i<sz; ++i) r[i] = a[i];
}

  //  a -> r: double to int16 
inline void 
  ZpzDom<Std16>::d2i ( const size_t sz, Array r, const double* a ) const
{
  union d_2_l {
    double d;
    int32 r[2];
  };
//  static const double offset = 4503599627370496.0; // 2^52
  double offset = 4503599627370496.0; // 2^52
  for (size_t i=0; i<sz; ++i)
  {
      register d_2_l tmp;
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
inline std::istream& ZpzDom<Std16>::read (std::istream& s) 
{
  char ch; 
  s >> std::ws >> ch;
  if (ch != '(')
//    GivError::throw_error( GivBadFormat("ZpzDom<Std16>::read: syntax error: no '('"));
    std::cerr << "GivBadFormat(ZpzDom<Std16>::read: syntax error: no '('))" << std::endl;

  s >> std::ws >> ch;
  if (ch != 'z')
//    GivError::throw_error( GivBadFormat("ZpzDom<Std16>::read: bad domain object"));
    std::cerr << "GivBadFormat(ZpzDom<Std16>::read: bad domain object))" << std::endl;

  s >> std::ws >> ch;
  if (ch != ',')
//    GivError::throw_error( GivBadFormat("ZpzDom<Std16>::read: syntax error: no ','"));
    std::cerr << "GivBadFormat(ZpzDom<Std16>::read: syntax error: no ',')) " << std::endl;

  s >> std::ws >> _p;

  s >> std::ws >> ch;
  if (ch != ')')
//    GivError::throw_error( GivBadFormat("ZpzDom<Std16>::read: syntax error: no ')'"));
    std::cerr << "GivBadFormat(ZpzDom<Std16>::read: syntax error: no ')')) " << std::endl;

  return s;
}

inline std::ostream& ZpzDom<Std16>::write (std::ostream& s ) const
{
  return s << "Std16 Givaro Z/pZ modulo " << residu();
}

inline std::istream& ZpzDom<Std16>::read (std::istream& s, Rep& a) const
{
  s >> a;
  init(a, a);
  return s;
}

inline std::ostream& ZpzDom<Std16>::write (std::ostream& s, const Rep a) const
{
  return s << a;
}
