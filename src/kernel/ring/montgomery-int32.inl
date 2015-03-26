// ==========================================================================
// $Id: givmontg32.inl,v 1.14 2011-02-04 14:11:46 jgdumas Exp $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// ==========================================================================

#ifndef __GIVARO_mong32_INL
#define __GIVARO_mong32_INL

namespace Givaro {

inline Montgomery<int32_t>::Element Montgomery<int32_t>::redcal(const Element c) const
{
        Element c0 = (Element)(c & MASK32);	/* c mod B */
        c0 = (Element)((c0 * _nim) & MASK32); 	/* -c/p mod B */
	// c0 *= _p;
	// c0 += c;			/* c = 0 mod B */
	c0 = c + c0*_p;
	c0 >>= HALF_BITS32;
        return (c0>_p?c0-=_p:c0);
}
inline Montgomery<int32_t>::Element Montgomery<int32_t>::redcsal(const Element c) const
{
        Element c0 = (Element)((c * _nim) & MASK32); 	/* -c/p mod B */
        c0 = c + c0 * _p; 		/* c = 0 mod B */
        c0 >>= HALF_BITS32;
	return (c0>_p?c0-=_p:c0);
}

inline Montgomery<int32_t>::Element& Montgomery<int32_t>::redc(Element& r, const Element c) const
{
        r = (Element)(c & MASK32);			/* c mod B */
	r *= _nim;
	r &= MASK32;
	r *= _p;
	r += c;
        r >>= HALF_BITS32;
	return (r>_p?r-=_p:r);
}

inline Montgomery<int32_t>::Element& Montgomery<int32_t>::redcs(Element& r, const Element c) const
{
        r = (Element)((c * _nim) & MASK32); 	/* -c/p mod B */
        r = c + r * _p; 		/* c = 0 mod B */
        r >>= HALF_BITS32;
	return (r>_p?r-=_p:r);
}

inline Montgomery<int32_t>::Element& Montgomery<int32_t>::redcin(Element& r) const
{
        Element c0 = (Element)(r & MASK32);	/* c mod B */
        c0 = (Element)((c0 * _nim) & MASK32); 	/* -c/p mod B */
        r += c0 * _p; 			/* c = 0 mod B */
        r >>= HALF_BITS32;
	return (r>_p?r-=_p:r);
}
inline Montgomery<int32_t>::Element& Montgomery<int32_t>::redcsin(Element& r) const
{
        Element c0 = (Element)((r * _nim) & MASK32); 	/* -c/p mod B */
        r += c0 * _p; 				/* c = 0 mod B */
        r >>= HALF_BITS32;
	return (r>_p?r-=_p:r);
}

} // namespace Givaro

// r = a*b
#define __GIVARO_MONTG32_MUL(r,p,a,b) (redc(r,a*b))
// r *= a
#define __GIVARO_MONTG32_MULIN(r,p,a) (redcin(r*=a))

// r = a - b
#define __GIVARO_MONTG32_SUB(r,p,a,b) ( r = (a>=b)? a-b: (p-b)+a )
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

#define __GIVARO_MONTG32_NEG(r,p,a) (r = (a == 0 ? 0 : p-a))
#define __GIVARO_MONTG32_NEGIN(r,p) (r = (r == 0 ? 0 : p-r))

namespace Givaro {

inline Montgomery<int32_t>::Residu_t Montgomery<int32_t>::residu( ) const
{ return _p; }

inline Montgomery<int32_t>::Rep& Montgomery<int32_t>::mul (Rep& r, const Rep a, const Rep b) const
{
  return __GIVARO_MONTG32_MUL(r,_p,a,b);
}

inline Montgomery<int32_t>::Rep& Montgomery<int32_t>::sub (Rep& r, const Rep a, const Rep b) const
{
  return __GIVARO_MONTG32_SUB(r,_p,a,b);
}

inline Montgomery<int32_t>::Rep& Montgomery<int32_t>::add (Rep& r, const Rep a, const Rep b) const
{
  __GIVARO_MONTG32_ADD(r,_p,a,b);
  return r;
}

inline Montgomery<int32_t>::Rep& Montgomery<int32_t>::neg (Rep& r, const Rep a) const
{
  return __GIVARO_MONTG32_NEG(r,_p,a);
}

inline Montgomery<int32_t>::Rep& Montgomery<int32_t>::inv (Rep& r, const Rep a) const
{
	// invext(aB) --> 1/a*1/B
	// % * B^3    --> B²/a
	// redc       --> B/a
    int32_t t;
    invext(t, int32_t(a), int32_t(_p));
    if (t < 0) t += _p;
    return redc(r, uint32_t(t) * _B3p) ;
}

inline Montgomery<int32_t>::Rep& Montgomery<int32_t>::div (Rep& r, const Rep a, const Rep b) const
{
	return mulin( inv(r,b), a );
}



inline Montgomery<int32_t>::Rep& Montgomery<int32_t>::mulin (Rep& r, const Rep a) const
{
  return __GIVARO_MONTG32_MULIN(r,_p, a);
}

inline Montgomery<int32_t>::Rep& Montgomery<int32_t>::divin (Rep& r, const Rep a) const
{
  Montgomery<int32_t>::Rep ia;
  inv(ia, a);
  return mulin(r, ia);
}

inline Montgomery<int32_t>::Rep& Montgomery<int32_t>::addin (Rep& r, const Rep a) const
{
   __GIVARO_MONTG32_ADDIN(r,_p, a);
  return r;
}

inline Montgomery<int32_t>::Rep& Montgomery<int32_t>::subin (Rep& r, const Rep a) const
{
  __GIVARO_MONTG32_SUBIN(r,_p, a);
  return r;
}


inline Montgomery<int32_t>::Rep& Montgomery<int32_t>::negin (Rep& r) const
{
  return __GIVARO_MONTG32_NEGIN(r,_p);
}

inline Montgomery<int32_t>::Rep& Montgomery<int32_t>::invin (Rep& r) const
{
	uint32_t t;
	return r = inv(t,r);
}

inline Montgomery<int32_t>::Rep& Montgomery<int32_t>::axpy
 (Rep& r, const Rep a, const Rep b, const Rep c) const
{
  __GIVARO_MONTG32_MULADD(r, _p, a, b, c);
  return r;
}

inline Montgomery<int32_t>::Rep&  Montgomery<int32_t>::axpyin
 (Rep& r, const Rep a, const Rep b) const
{
  __GIVARO_MONTG32_MULADDIN(r, _p, a, b);
  return r;
}

// r <- a*b-c
inline Montgomery<int32_t>::Rep&  Montgomery<int32_t>::axmy
 (Rep& r, const Rep a, const Rep b, const Rep c) const
{
    return this->subin(this->mul(r,a,b), c);
}

// r = c - a*b
inline Montgomery<int32_t>::Rep&  Montgomery<int32_t>::maxpy
 (Rep& r, const Rep a, const Rep b, const Rep c) const
{
    Rep t;
    return this->sub(r, c, this->mul(t,a,b));
}

// r -= a*b
inline Montgomery<int32_t>::Rep&  Montgomery<int32_t>::maxpyin
 (Rep& r, const Rep a, const Rep b) const
{
    Rep t;
    return this->subin(r, this->mul(t,a,b));
}

// r = a*b - r
inline Montgomery<int32_t>::Rep&  Montgomery<int32_t>::axmyin
 (Rep& r, const Rep a, const Rep b) const
{
    maxpyin(r,a,b);
    return negin(r);
}



 // ------------------------- Miscellaneous functions

inline int Montgomery<int32_t>::isZero(const Rep a) const
{ return a == Montgomery<int32_t>::zero; }

inline int Montgomery<int32_t>::isOne(const Rep a) const
{ return a == Montgomery<int32_t>::one; }

inline int Montgomery<int32_t>::isMOne(const Rep a) const
{ return a == Montgomery<int32_t>::mOne; }



inline size_t Montgomery<int32_t>::length(const Rep) const
{ return Montgomery<int32_t>::size_rep;}

// ---------
// -- misc operations
// ---------


inline  Montgomery<int32_t>::Rep&  Montgomery<int32_t>::init ( Rep& r, const double a ) const
{
  int sign; double ua;
  if (a < 0.0) { sign =-1; ua = -a;}
  else { ua = a; sign =1; }
  if ( ua > Signed_Trait<uint32_t>::max()){
//    ua -= (double)floor(ua * _invdp)*_dp;
    ua = fmod(ua,_dp);
    r = (Rep) ua;
  } else
    r = (ua >=_p) ? (uint32_t) ua % (uint32_t)_p : (uint32_t) ua;
  if (r && (sign ==-1)) r = _p - r;
//  std::cerr << a << "dbl --> " << r << "(" << redcal(r*_B2p) << ")" << std::endl;
  return redc(r,r*_B2p);
}

inline  Montgomery<int32_t>::Rep&  Montgomery<int32_t>::init ( Rep& r, const float a ) const
{
    return init(r, (double)a);
}



inline  Montgomery<int32_t>::Rep&  Montgomery<int32_t>::init ( Rep& r, const unsigned long a ) const
{
	r = Rep( a >= (uint32_t)_p ? a % (uint32_t)_p : a);
  return redc(r,r*_B2p);
}

inline  Montgomery<int32_t>::Rep&  Montgomery<int32_t>::init ( Rep& r, const long a ) const
{
  int sign; unsigned long ua;
  if (a <0)
  {
	  sign =-1;
	  ua = (unsigned long) -a;
  }
  else {
	  ua = (unsigned long)a;
	  sign =1;
  }
  r =Rep ( ua >= (uint32_t)_p ? ua % (uint32_t)_p : ua);
  if (r && (sign ==-1))
	  r = _p - r;
  return redc(r,r*_B2p);
}


inline  Montgomery<int32_t>::Rep&  Montgomery<int32_t>::init ( Rep& r, const Integer& Residu ) const
{
  long tr;
  if (Residu <0) {
      // -a = b [p]
      // a = p-b [p]
    if ( Residu <= (Integer)(-_p) ) tr = long( (-Residu) % _p) ;
    else tr = long(-Residu);
    if (tr)
      r = Rep(_p - (unsigned long)tr);
    else
      r = zero;
  } else {
    if (Residu >= (Integer)_p ) tr =   long(Residu % _p) ;
    else tr = long(Residu);
    r = Rep(tr);
  }
  return redc(r,r*_B2p);
}





inline Montgomery<int32_t>::Rep& Montgomery<int32_t>::init( Rep& a, const int i) const
{ return init(a,(long)i); }

inline Montgomery<int32_t>::Rep& Montgomery<int32_t>::init( Rep& a, const unsigned int i) const
{ return init(a,(unsigned long)i); }

inline  Montgomery<int32_t>::Rep&  Montgomery<int32_t>::assign
  ( Rep& r, const Rep a ) const
{ return r=a; }



inline Montgomery<int32_t>::Rep& Montgomery<int32_t>::init ( Rep& r ) const
{ return r = zero; }


template< class Random >
inline  Montgomery<int32_t>::Rep& Montgomery<int32_t>::random(Random& g, Rep& a) const
{
	        return init(a, g());
}

template< class Random >
inline  Montgomery<int32_t>::Rep& Montgomery<int32_t>::random(Random& g, Rep& a, const Rep& b) const
{
	        return init(a, g());
}
template< class Random >
inline  Montgomery<int32_t>::Rep& Montgomery<int32_t>::random(Random& g, Rep& a, long b) const
{
	        return init(a, g() %(uint32_t) b);

}

template< class Random >
inline  Montgomery<int32_t>::Rep& Montgomery<int32_t>::nonzerorandom(Random& g, Rep& a) const
{
	        while (isZero(init(a, g()))) {};
		return a;
}

template< class Random >
inline  Montgomery<int32_t>::Rep& Montgomery<int32_t>::nonzerorandom(Random& g, Rep& a, const Rep& b) const
{
	        while (isZero(init(a, g()))) {};
		return a;
}

template< class Random >
inline  Montgomery<int32_t>::Rep& Montgomery<int32_t>::nonzerorandom(Random& g, Rep& a, long b) const
{
	        while (isZero(init(a, g() %(uint32_t) b))) {};
		return a;
}


 // -- Input: (z, <_p>)
inline std::istream& Montgomery<int32_t>::read (std::istream& s)
{
  char ch;
  s >> std::ws >> ch;
  if (ch != '(')
//    GivError::throw_error( GivBadFormat("Montgomery<int32_t>::read: syntax error: no '('"));
    std::cerr << "GivBadFormat(Montgomery<int32_t>::read: syntax error: no '('))" << std::endl;

  s >> std::ws >> ch;
  if (ch != 'z')
//    GivError::throw_error( GivBadFormat("Montgomery<int32_t>::read: bad domain object"));
    std::cerr << "GivBadFormat(Montgomery<int32_t>::read: bad domain object))" << std::endl;

  s >> std::ws >> ch;
  if (ch != ',')
//    GivError::throw_error( GivBadFormat("Montgomery<int32_t>::read: syntax error: no ','"));
    std::cerr << "GivBadFormat(Montgomery<int32_t>::read: syntax error: no ',')) " << std::endl;

  s >> std::ws >> _p;

  s >> std::ws >> ch;
  if (ch != ')')
//    GivError::throw_error( GivBadFormat("Montgomery<int32_t>::read: syntax error: no ')'"));
    std::cerr << "GivBadFormat(Montgomery<int32_t>::read: syntax error: no ')')) " << std::endl;

  return s;
}

inline std::ostream& Montgomery<int32_t>::write (std::ostream& s ) const
{
  return s << "Givaro Montgomery Z/pZ, p=" << residu();
}

inline std::istream& Montgomery<int32_t>::read (std::istream& s, Rep& a) const
{
  Integer tmp;
  s >> tmp;
  init(a, tmp);
  return s;
}

inline std::ostream& Montgomery<int32_t>::write (std::ostream& s, const Rep a) const
{
    Rep tmp;
    return s << redcs(tmp,a);
}

} // namespace Givaro

#endif // __GIVARO_mong32_INL
