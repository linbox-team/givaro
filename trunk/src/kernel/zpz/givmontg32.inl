// ==========================================================================
// $Id: givmontg32.inl,v 1.4 2004-10-11 12:29:50 jgdumas Exp $
// ==========================================================================

inline Montgomery<Std32>::Element Montgomery<Std32>::redcal(const Element c) const {
        Element c0 = c & MASK32;	/* c mod B */
        c0 = (c0 * _nim) & MASK32; 	/* -c/p mod B */
	// c0 *= _p;
	// c0 += c;			/* c = 0 mod B */
	c0 = c + c0*_p;
	c0 >>= HALF_BITS32;
        return (c0>_p?c0-=_p:c0);
}
inline Montgomery<Std32>::Element Montgomery<Std32>::redcsal(const Element c) const {
        Element c0 = (c * _nim) & MASK32; 	/* -c/p mod B */
        c0 = c + c0 * _p; 		/* c = 0 mod B */
        c0 >>= HALF_BITS32;
	return (c0>_p?c0-=_p:c0);
}

inline Montgomery<Std32>::Element& Montgomery<Std32>::redc(Element& r, const Element c) const {
        r = c & MASK32;			/* c mod B */
	r *= _nim;
	r &= MASK32;
	r *= _p;
	r += c;
        r >>= HALF_BITS32;
	return (r>_p?r-=_p:r);
}

inline Montgomery<Std32>::Element& Montgomery<Std32>::redcs(Element& r, const Element c) const {
        r = (c * _nim) & MASK32; 	/* -c/p mod B */
        r = c + r * _p; 		/* c = 0 mod B */
        r >>= HALF_BITS32;
	return (r>_p?r-=_p:r);
}

inline Montgomery<Std32>::Element& Montgomery<Std32>::redcin(Element& r) const {
        Element c0 = r & MASK32;	/* c mod B */
        c0 = (c0 * _nim) & MASK32; 	/* -c/p mod B */
        r += c0 * _p; 			/* c = 0 mod B */
        r >>= HALF_BITS32;
	return (r>_p?r-=_p:r);
}
inline Montgomery<Std32>::Element& Montgomery<Std32>::redcsin(Element& r) const {
        Element c0 = (r * _nim) & MASK32; 	/* -c/p mod B */
        r += c0 * _p; 				/* c = 0 mod B */
        r >>= HALF_BITS32;
	return (r>_p?r-=_p:r);
}

// r = a*b
#define __GIVARO_MONTG32_MUL(r,p,a,b) (redc(r,a*b))
// r *= a
#define __GIVARO_MONTG32_MULIN(r,p,a) (redcin(r*=a))

// r = a - b
#define __GIVARO_MONTG32_SUB(r,p,a,b) ( r = (a>b)? a-b: (p-b)+a )
// r -= a
//#define __GIVARO_MONTG32_SUBIN(r,p,a) { r -= a; r= (r < 0 ? r+p : r); }
#define __GIVARO_MONTG32_SUBIN(r,p,a) { if (r<a) r+=(p-a); else r-=a; }

// r = a+b
#define __GIVARO_MONTG32_ADD(r,p,a,b) { r = (a+b); r= (r < p ? r : r-p); }
// r += a
#define __GIVARO_MONTG32_ADDIN(r,p,a) { r += a;  r= (r < p ? r : r-p); }

// r <- a*b+c % p
//#define __GIVARO_MONTG32_MULADD(r,p,a,b,c) { redc(r,a*b) += c; r= (r < p ? r : r-p); }
#define __GIVARO_MONTG32_MULADD(r,p,a,b,c) { r = redcal(a*b) + c; r= (r < p ? r : r-p); }


#define __GIVARO_MONTG32_MULADDIN(r,p,a,b) { r+=redcal(a*b); r= (r < p ? r : r-p); }

// a*b-c
#define __GIVARO_MONTG32_MULSUB(r,p,a,b,c) (r)
// a*b-c
#define __GIVARO_MONTG32_SUBMULIN(r,p,a,b) (r)

#define __GIVARO_MONTG32_NEG(r,p,a) (r = (a == 0 ? 0 : p-a))
#define __GIVARO_MONTG32_NEGIN(r,p) (r = (r == 0 ? 0 : p-r))


inline Montgomery<Std32>::Residu_t Montgomery<Std32>::residu( ) const
{ return _p; }

inline Montgomery<Std32>::Rep& Montgomery<Std32>::mul (Rep& r, const Rep a, const Rep b) const
{ 
  return __GIVARO_MONTG32_MUL(r,_p,a,b); 
}

inline Montgomery<Std32>::Rep& Montgomery<Std32>::sub (Rep& r, const Rep a, const Rep b) const
{ 
  return __GIVARO_MONTG32_SUB(r,_p,a,b); 
}

inline Montgomery<Std32>::Rep& Montgomery<Std32>::add (Rep& r, const Rep a, const Rep b) const
{ 
  __GIVARO_MONTG32_ADD(r,_p,a,b); 
  return r; 
}

inline Montgomery<Std32>::Rep& Montgomery<Std32>::neg (Rep& r, const Rep a) const
{ 
  return __GIVARO_MONTG32_NEG(r,_p,a); 
}

inline Montgomery<Std32>::Rep& Montgomery<Std32>::inv (Rep& r, const Rep a) const
{ 
	// invext(aB) --> 1/a*1/B
	// % * B^3    --> B²/a
	// redc       --> B/a
	int32 t;
	return redc(r, uint32( invext( t,int32(a),int32(_p)) ) * _B3p) ;
}

inline Montgomery<Std32>::Rep& Montgomery<Std32>::div (Rep& r, const Rep a, const Rep b) const
{ 
	return mulin( inv(r,b), a );
}



inline Montgomery<Std32>::Rep& Montgomery<Std32>::mulin (Rep& r, const Rep a) const
{ 
  return __GIVARO_MONTG32_MULIN(r,_p, a); 
}

inline Montgomery<Std32>::Rep& Montgomery<Std32>::divin (Rep& r, const Rep a) const
{ 
  Montgomery<Std32>::Rep ia;
  inv(ia, a);
  return mulin(r, ia);
}

inline Montgomery<Std32>::Rep& Montgomery<Std32>::addin (Rep& r, const Rep a) const
{ 
   __GIVARO_MONTG32_ADDIN(r,_p, a); 
  return r; 
}

inline Montgomery<Std32>::Rep& Montgomery<Std32>::subin (Rep& r, const Rep a) const
{ 
  __GIVARO_MONTG32_SUBIN(r,_p, a); 
  return r; 
}


inline Montgomery<Std32>::Rep& Montgomery<Std32>::negin (Rep& r) const
{ 
  return __GIVARO_MONTG32_NEGIN(r,_p); 
}

inline Montgomery<Std32>::Rep& Montgomery<Std32>::invin (Rep& r) const
{ 
	uint32 t;
	return r = inv(t,r);
}

inline Montgomery<Std32>::Rep& Montgomery<Std32>::axpy 
 (Rep& r, const Rep a, const Rep b, const Rep c) const
{ 
  __GIVARO_MONTG32_MULADD(r, _p, a, b, c); 
  return r;
}

inline Montgomery<Std32>::Rep&  Montgomery<Std32>::axpyin 
 (Rep& r, const Rep a, const Rep b) const
{ 
  __GIVARO_MONTG32_MULADDIN(r, _p, a, b); 
  return r;
}

inline Montgomery<Std32>::Rep&  Montgomery<Std32>::axmy
 (Rep& r, const Rep a, const Rep b, const Rep c) const
{
  return __GIVARO_MONTG32_MULSUB(r, _p, a, b, c);
}

// r -= a*b
inline Montgomery<Std32>::Rep&  Montgomery<Std32>::axmyin 
 (Rep& r, const Rep a, const Rep b) const
{
  return __GIVARO_MONTG32_SUBMULIN(r, _p, a, b );
}



 // ------------------------- Miscellaneous functions

inline int Montgomery<Std32>::iszero(const Rep a) const
{ return a == Montgomery<Std32>::zero; }

inline int Montgomery<Std32>::isone(const Rep a) const
{ return a == Montgomery<Std32>::one; }

inline int Montgomery<Std32>::isZero( const Rep a ) const {
  return iszero(a);
}

inline int Montgomery<Std32>::isOne ( const Rep a ) const {
//	std::cerr << "one: " << one << ", a: " << a << std::endl;
  return isone(a);
}


inline size_t Montgomery<Std32>::length(const Rep a) const
{ return Montgomery<Std32>::size_rep;}

// ---------
// -- misc operations
// ---------


inline  Montgomery<Std32>::Rep&  Montgomery<Std32>::init ( Rep& r, const double a ) const
{
  int sign; double ua;
  if (a < 0.0) { sign =-1; ua = -a;}
  else { ua = a; sign =1; }
  if ( ua > Signed_Trait<uint32>::max()){
    ua -= (double)floor(ua * _invdp)*_dp; 
    r = (Rep) ua;
  } else
    r = (ua >=_p) ? (uint32) ua % (uint32)_p : (uint32) ua;
  if (r && (sign ==-1)) r = _p - r;
//  std::cerr << a << "dbl --> " << r << "(" << redcal(r*_B2p) << ")" << std::endl;
  return redc(r,r*_B2p);
}

inline  Montgomery<Std32>::Rep&  Montgomery<Std32>::init ( Rep& r, const float a ) const {
    return init(r, (double)a);
}



inline  Montgomery<Std32>::Rep&  Montgomery<Std32>::init ( Rep& r, const unsigned long a ) const
{ r = ( a >= (uint32)_p ? a % (uint32)_p : a);
  return redc(r,r*_B2p);
}

inline  Montgomery<Std32>::Rep&  Montgomery<Std32>::init ( Rep& r, const long a ) const
{
  int sign; unsigned long ua;
  if (a <0) { sign =-1; ua = -a;}
  else { ua = a; sign =1; }
  r = ( ua >= (uint32)_p ? ua % (uint32)_p : ua);
  if (r && (sign ==-1)) r = _p - r;
  return redc(r,r*_B2p);
}


inline  Montgomery<Std32>::Rep&  Montgomery<Std32>::init ( Rep& r, const Integer& residu ) const 
{
  long tr;
  if (residu <0) {
      // -a = b [p]
      // a = p-b [p]
    if ( residu <= (Integer)(-_p) ) tr = long( (-residu) % _p) ;
    else tr = long(-residu);
    if (tr)
      r = _p - (unsigned long)tr;
    else
      r = zero;
  } else {
    if (residu >= (Integer)_p ) tr =   long(residu % _p) ;
    else tr = long(residu);
    r = tr;
  }
  return redc(r,r*_B2p);
}





inline Montgomery<Std32>::Rep& Montgomery<Std32>::init( Rep& a, const int i) const { return init(a,(long)i); }

inline Montgomery<Std32>::Rep& Montgomery<Std32>::init( Rep& a, const unsigned int i) const { return init(a,(unsigned long)i); }

inline  Montgomery<Std32>::Rep&  Montgomery<Std32>::assign 
  ( Rep& r, const Rep a ) const
{ return r=a; }



inline Montgomery<Std32>::Rep& Montgomery<Std32>::init ( Rep& r ) const
{ return r = zero; }


template< class RandIter >
inline  Montgomery<Std32>::Rep& Montgomery<Std32>::random(RandIter& g, Rep& a) const {
	        return init(a, g());
}

template< class RandIter >
inline  Montgomery<Std32>::Rep& Montgomery<Std32>::random(RandIter& g, Rep& a, const Rep& b) const {
	        return init(a, g());
}
template< class RandIter >
inline  Montgomery<Std32>::Rep& Montgomery<Std32>::random(RandIter& g, Rep& a, long b) const {
	        return init(a, g() %(uint32) b);
		
}

template< class RandIter >
inline  Montgomery<Std32>::Rep& Montgomery<Std32>::nonzerorandom(RandIter& g, Rep& a) const {
	        while (iszero(init(a, g()))) {};
		return a;
}

template< class RandIter >
inline  Montgomery<Std32>::Rep& Montgomery<Std32>::nonzerorandom(RandIter& g, Rep& a, const Rep& b) const {
	        while (iszero(init(a, g()))) {};
		return a;
}

template< class RandIter >
inline  Montgomery<Std32>::Rep& Montgomery<Std32>::nonzerorandom(RandIter& g, Rep& a, long b) const {
	        while (iszero(init(a, g() %(uint32) b))) {};
		return a;
}


 // -- Input: (z, <_p>)
inline std::istream& Montgomery<Std32>::read (std::istream& s) 
{
  char ch; 
  s >> std::ws >> ch;
  if (ch != '(')
//    GivError::throw_error( GivBadFormat("Montgomery<Std32>::read: syntax error: no '('"));
    std::cerr << "GivBadFormat(Montgomery<Std32>::read: syntax error: no '('))" << std::endl;

  s >> std::ws >> ch;
  if (ch != 'z')
//    GivError::throw_error( GivBadFormat("Montgomery<Std32>::read: bad domain object"));
    std::cerr << "GivBadFormat(Montgomery<Std32>::read: bad domain object))" << std::endl;

  s >> std::ws >> ch;
  if (ch != ',')
//    GivError::throw_error( GivBadFormat("Montgomery<Std32>::read: syntax error: no ','"));
    std::cerr << "GivBadFormat(Montgomery<Std32>::read: syntax error: no ',')) " << std::endl;

  s >> std::ws >> _p;

  s >> std::ws >> ch;
  if (ch != ')')
//    GivError::throw_error( GivBadFormat("Montgomery<Std32>::read: syntax error: no ')'"));
    std::cerr << "GivBadFormat(Montgomery<Std32>::read: syntax error: no ')')) " << std::endl;

  return s;
}

inline std::ostream& Montgomery<Std32>::write (std::ostream& s ) const
{
  return s << "Givaro Montgomery Z/pZ, p=" << residu();
}

inline std::istream& Montgomery<Std32>::read (std::istream& s, Rep& a) const
{
  s >> a;
  init(a, a);
  return s;
}

inline std::ostream& Montgomery<Std32>::write (std::ostream& s, const Rep a) const
{
    Rep tmp;
    return s << redcs(tmp,a);
}
