// ========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int.inl,v $
// Copyright(c)'1994-2010 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: M. Samama, T. Gautier
// Modified by: B. Boyer
// $Id: gmp++_int.inl,v 1.19 2011-01-19 18:29:09 bboyer Exp $
// ========================================================================
// Description:

#define GMP__ABS(l)     ((l) <0 ? -l : l)
#define GMP__SGN(l)    ((l) <0 ? -1 : (l >0 ? 1 : 0))

#include <cassert>
#include <givaro/givtimer.h>

//-----------------------------~Integer()
inline Integer::~Integer() {  mpz_clear((mpz_ptr)&gmp_rep) ; }

//-------------------------------Integer(const Integer &n)
inline Integer::Integer(const Integer &n) {
	mpz_init_set ( (mpz_ptr)&gmp_rep, (mpz_ptr)&(n.gmp_rep)) ;
}

//------------------------------------------operator = (const Integer &n)
inline Integer& Integer::logcpy(const Integer &n)
{
	if (this == &n) return *this;
	mpz_set ( (mpz_ptr)&gmp_rep, (mpz_ptr)&(n.gmp_rep)) ;
	return *this;
}

// same as logcopy
inline Integer& Integer::operator = (const Integer &n) { return logcpy(n) ; }

//-----------------------------Integer(int n)
inline Integer::Integer(int n) { mpz_init_set_si((mpz_ptr)&gmp_rep, n) ; }

//-----------------------------Integer(uint n)
inline Integer::Integer(unsigned char n) { mpz_init_set_ui((mpz_ptr)&gmp_rep, n) ; }

//-----------------------------Integer(uint n)
inline Integer::Integer(unsigned int n) { mpz_init_set_ui((mpz_ptr)&gmp_rep, n) ; }

//-----------------------------Integer(long n)
inline Integer::Integer(long n) { mpz_init_set_si((mpz_ptr)&gmp_rep, n) ; }

//-----------------------------Integer(unsigned long n)
inline Integer::Integer(unsigned long n) { mpz_init_set_ui((mpz_ptr)&gmp_rep, n) ; }

#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
#include <stdio.h>
//-----------------------------Integer(long long n)
// log[10](2^8) < 2.408239966
inline Integer::Integer(long long n) {
	char * tmp = new char[long(2.408239966*sizeof(long long))+1]; sprintf(tmp,"%lld",n);
	mpz_init_set_str((mpz_ptr)&gmp_rep, tmp, 10) ;
	delete [] tmp;
}

//-----------------------------Integer(unsigned long long n)
// log[10](2^8) < 2.408239966
inline Integer::Integer(unsigned long long n) {
	char * tmp = new char[ long(2.408239966*sizeof(unsigned long long))+1];
	sprintf(tmp,"%llu",n);
	mpz_init_set_str((mpz_ptr)&gmp_rep, tmp, 10) ;
	delete [] tmp;
}
#endif


//-----------------------------Integer(double)
inline Integer::Integer(double d) { mpz_init_set_d((mpz_ptr)&gmp_rep, d) ; }


//-----------------------------Integer(const neutral n), default n = zero
/* Neutral is causing a problem
   inline Integer::Integer(const Neutral n) {
   if (n == Neutral::zero) mpz_init_set_ui((mpz_ptr)&gmp_rep, 0L) ;
   else  mpz_init_set_ui((mpz_ptr)&gmp_rep, 1L) ;
   }
   */

//-------------------------------------------------inline comparaison operators
inline int operator != (const Integer& a , const Integer& b)
{ return compare(a,b) != 0; }

inline int operator != (int l, const Integer& n)
{ return n.operator != (l); }

inline int operator != (long l, const Integer& n)
{ return n.operator != (l); }

inline int operator != (unsigned long l, const Integer& n)
{ return n.operator != (l); }

inline int operator == (const Integer& a, const Integer& b)
{  return compare(a,b) == 0; }

inline int operator == (int l, const Integer& n)
{ return (! (n.operator != (l))); }

inline int operator == (long l, const Integer& n)
{ return (! (n.operator != (l))); }

inline int operator == (unsigned long l, const Integer& n)
{ return (! (n.operator != (l))); }

inline int operator == (const Integer& n, unsigned long l)
{ return (! (n.operator != (l))); }

inline int operator == (const Integer& n, int l)
{ return (! (n.operator != (l))); }

inline int operator == (const Integer& n, long l)
{ return (! (n.operator != (l))); }

inline int operator < (const Integer& a , const Integer& b)
{ return compare(a,b) < 0; }

inline int operator < (const int l, const Integer& n)
{ return n > l; }

inline int operator < (const long l, const Integer& n)
{ return n > l; }

inline int operator < (const unsigned long l, const Integer& n)
{ return n > l; }

inline int operator <= (const Integer& n, unsigned long l)
{  return (! (n > l) ); }

inline int operator <= (unsigned long l, const Integer& n)
{  return (! (n < l) );}

inline int operator >= (unsigned long l, const Integer& n)
{  return (! (n < l) );}

inline int operator >= (const Integer& n, unsigned long l)
{  return (! (n < l) );}

inline int operator > (int l, const Integer& n)
{ return n < l; }

inline int operator > (long l, const Integer& n)
{ return n < l; }

inline int operator > (unsigned long l, const Integer& n)
{ return n < l; }

inline int operator >  (const Integer& a , const Integer& b)
{ return compare(a,b) > 0; }

inline int operator <= (const Integer& a, const Integer& b)
{ return compare(a,b) <= 0; }

inline int operator <= (const Integer& n, int l)
{  return (! (n > l) ); }

inline int operator <= (const Integer& n, long l)
{  return (! (n > l) ); }

inline int operator <= (int l, const Integer& n)
{  return (! (n < l) );}

inline int operator <= (long l, const Integer& n)
{  return (! (n < l) );}

inline int operator >= (const Integer& a, const Integer& b)
{ return compare(a,b) >= 0; }

inline int operator >= (int l, const Integer& n)
{  return (! (n > l) );}

inline int operator >= (long l, const Integer& n)
{  return (! (n > l) );}

inline int operator >= (const Integer& n, int l)
{  return (! (n < l) );}

inline int operator >= (const Integer& n, long l)
{  return (! (n < l) );}


//----------------------------------arithmetic inline operators
inline Integer Integer::operator - () const
{
	// JGD 18.06.1999
	Integer Res ;
	mpz_neg((mpz_ptr)&Res.gmp_rep, (mpz_ptr)&gmp_rep );
	return Res ;
}

// -- operator +
inline Integer operator + (const int l, const Integer& n) { return n + (long)l; }
inline Integer operator + (const unsigned int l, const Integer& n) { return n + (unsigned long)l; }
inline Integer operator + (const long l, const Integer& n) { return n + l; }
inline Integer operator + (const unsigned long l, const Integer& n) { return n + l; }
inline Integer operator + (const Integer& n, const int l) { return n + (long)l; }
inline Integer operator + (const Integer& n, const unsigned int l) { return n + (unsigned long)l; }

inline Integer& operator += (Integer& n, const int l) { return n += (long)l; }
inline Integer& operator += (Integer& n, const unsigned int l) { return n += (unsigned long)l; }

#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
inline Integer operator + (const Integer& n, const long long l) {return n + (Integer)l; }
inline Integer operator + (const Integer& n, const unsigned long long l) {return n + (Integer)l; }
inline Integer operator + (const long long l, const Integer& n) {return n+l;}
inline Integer operator + (const unsigned long long l, const Integer& n) {return n+l;}
inline Integer& operator += (Integer& n, const long long l) { return n += (Integer)l; }
inline Integer& operator += (Integer& n, const unsigned long long l) { return n += (Integer)l; }
#endif


// -- operator -
inline Integer operator - (const int l, const Integer& n) { return -(n - (long)l); }
inline Integer operator - (const unsigned int l, const Integer& n) { return -(n - (unsigned long)l); }
inline Integer operator - (const long l, const Integer& n) { return -(n - l); }
inline Integer operator - (const unsigned long l, const Integer& n) { return -(n - l); }
inline Integer operator - (const Integer& n, const int l) { return n - (long)l; }
inline Integer operator - (const Integer& n, const unsigned int l) { return n - (unsigned long)l; }

inline Integer& operator -= (Integer& n, const int l) { return n -= (long)l; }
inline Integer& operator -= (Integer& n, const unsigned int l) { return n -= (unsigned long)l; }

#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
inline Integer operator - (const Integer& n, const long long l) {return n - (Integer)l; }
inline Integer operator - (const Integer& n, const unsigned long long l) {return n - (Integer)l; }
inline Integer operator - (const long long l, const Integer& n) {return n-l;}
inline Integer operator - (const unsigned long long l, const Integer& n) {return n-l;}
inline Integer& operator -= (Integer& n, const long long l) { return n -= (Integer)l; }
inline Integer& operator -= (Integer& n, const unsigned long long l) { return n -= (Integer)l; }
#endif

// -- operator *
inline Integer operator * (const int l, const Integer& n) { return n * (long)l; }
inline Integer operator * (const unsigned int l, const Integer& n) { return n * (unsigned long)l; }
inline Integer operator * (const long l, const Integer& n) { return n * l; }
inline Integer operator * (const unsigned long l, const Integer& n) { return n * l; }
inline Integer operator * (const Integer& n, const int l) { return n * (long)l; }
inline Integer operator * (const Integer& n, const unsigned int l) { return n * (unsigned long)l; }

inline Integer& operator *= (Integer& n, const int l) { return n *= (long)l; }
inline Integer& operator *= (Integer& n, const unsigned int l) { return n *= (unsigned long)l; }

#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
inline Integer operator * (const Integer& n, const long long l) {return n * (Integer)l; }
inline Integer operator * (const Integer& n, const unsigned long long l) {return n * (Integer)l; }
inline Integer operator * (const long long l, const Integer& n) {return n*l;}
inline Integer operator * (const unsigned long long l, const Integer& n) {return n*l;}
inline Integer& operator *= (Integer& n, const long long l) { return n *= (Integer)l; }
inline Integer& operator *= (Integer& n, const unsigned long long l) { return n *= (Integer)l; }
#endif

// -- operator /
inline Integer operator / (const int l, const Integer& n) { return Integer(l)/n; }
inline Integer operator / (const long l, const Integer& n) { return Integer(l)/n; }
inline Integer operator / (const Integer& n, const int l) { return n / (long)l; }
inline Integer operator / (const Integer& n, const unsigned int l) { return n / (unsigned long)l; }

inline Integer& operator /= (Integer& n, const int l) { if (l>=0) return n /= (unsigned long)l; else return  n = -(n / (unsigned long)-l); }
inline Integer& operator /= (Integer& n, const long l) { return n /= (unsigned long)l; }
inline Integer& operator /= (Integer& n, const unsigned int l) { return n /= (unsigned long)l; }

// -- operator %
inline Integer operator % (const int l, const Integer& n) { return Integer(l) % n; }
inline Integer operator % (const long l, const Integer& n) { return Integer(l) % n; }
inline Integer operator % (const Integer& n, const int l) { return n % (long)l; }
inline Integer operator % (const Integer& n, const unsigned int l) { return n % (unsigned long)l; }

inline Integer& operator %= (Integer& n, const int l) { return n %= (long)l; }
inline Integer& operator %= (Integer& n, const unsigned int l) { return n %= (unsigned long)l; }


//----------miscellaneous inline functions

inline int Integer::priv_sign() const { return mpz_sgn( (mpz_ptr)&gmp_rep ); }

inline int isOne(const Integer& a) { return ! mpz_cmp_ui((mpz_ptr)&(a.gmp_rep), 1UL); }

inline int isZero(const Integer& a) { return ! mpz_cmp_ui((mpz_ptr)&(a.gmp_rep), 0UL); }

inline int isZero(const short int a) { return a ==0; }
inline int isZero(const int a) { return a ==0; }
inline int isZero(const long a) { return a ==0; }
inline int isZero(const unsigned short int a) { return a ==0; }
inline int isZero(const unsigned int a) { return a ==0; }
inline int isZero(const unsigned long a) { return a ==0UL; }
#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
inline int isZero(const unsigned long long a) { return a ==0ULL; }
inline int isZero(const long long a) { return a ==0LL; }
#endif

inline int sign(const Integer& a) { return a.priv_sign(); }

inline unsigned long length(const Integer& a) { return mpz_size( (mpz_ptr)&(a.gmp_rep) ) * sizeof(unsigned long); }

inline Integer abs(const Integer &n) { if (sign(n) >= 0) return n; return -n; }

inline size_t Integer::size() const { return  mpz_size( (mpz_ptr)&gmp_rep ) ; }

inline size_t Integer::size_in_base(int BASE) const { return  mpz_sizeinbase ((mpz_ptr)&gmp_rep, BASE);}

inline size_t Integer::bitsize() const { return  mpz_sizeinbase ((mpz_ptr)&gmp_rep, 2);}

inline unsigned long Integer::operator[](size_t i) const
{ if ( mpz_size( (mpz_ptr)&gmp_rep ) > i)
	return mpz_getlimbn( (mpz_ptr)&gmp_rep, i);
	else
		return 0;
}

//-------------------------------------------------inline >> & << operators
inline std::ostream& operator<< (std::ostream& o, const Integer& a) { return a.print(o); }

//-----------------------------------------------------
//----------------------- Random integers -------------
//-----------------------------------------------------

/* ********************** */
/* seeding, initialising  */
/* ********************** */
inline gmp_randclass& Integer::randstate()
{
	static gmp_randclass randstate(gmp_randinit_default);
	return static_cast<gmp_randclass&>(randstate);
}

inline void Integer::seeding(unsigned long  s)
{
	Integer::randstate().seed(s) ;
}

inline void Integer::seeding(Integer  s)
{
	Integer::randstate().seed((mpz_class) (mpz_ptr) &(s.gmp_rep) ) ;
}

inline void Integer::seeding()
{
	Integer::seeding( BaseTimer::seed() );
}


// BB : good seeding but not so efficient...
inline bool Integer::RandBool()
{
	if (Integer::random(1UL)) return true;
	else return false ;
}


/* ****************************** */
/*  random number smaller than m  */
/* ****************************** */

#ifdef __GMP_PLUSPLUS__
//! returns a random integer \p r in the intervall <code>[[0, m-1]]</code>
template<bool U>
inline Integer& Integer::random_lessthan (Integer& r, const Integer & m)
{
	mpz_set( (mpz_ptr) &(r.gmp_rep) ,
		 ( (mpz_class)Integer::randstate().get_z_range((mpz_class) (mpz_ptr) &(m.gmp_rep)) ).get_mpz_t() );
	if(!U) if (Integer::RandBool()) Integer::negin(r);
	return r;
}
#else
//! returns a random integer \p r in the intervall <code>[[0, m-1]]</code>
template<bool U>
inline Integer& Integer::random_lessthan (Integer& r, const Integer & m)
{
	mpz_urandomm((mpz_ptr) &(r.gmp_rep),Integer::randstate(),(mpz_ptr)&(m.gmp_rep));
	if(!U) if (Integer::RandBool()) Integer::negin(r);
	return r;
}
#endif

/* ******************************** */
/*  random number smaller than 2^m  */
/* ******************************** */

#ifdef __GMP_PLUSPLUS__
//! returns a random integer \p r of at most \p m bits
template<bool U>
inline Integer& Integer::random_lessthan_2exp (Integer& r, const unsigned long & m)
{
	mpz_set( (mpz_ptr) &(r.gmp_rep) , ((mpz_class)Integer::randstate().get_z_bits(m)).get_mpz_t() );
	if(!U) {if (Integer::RandBool()) Integer::negin(r);}
	return r;
}
#else
//! returns a random integer \p r of at most \p m bits
template<bool U>
inline Integer& Integer::random_lessthan_2exp (Integer& r, const unsigned long & m)
{
	mpz_urandomb((mpz_ptr) &(r.gmp_rep),Integer::randstate(),m) ;
	if(!U) if (Integer::RandBool()) Integer::negin(r);
	return r;
}
#endif

template<bool U>
inline Integer Integer::random_lessthan_2exp (const unsigned long & m)
{
       Integer r ;
       random_lessthan_2exp<U>(r,m);
       return r;
}

/* synonyms */
template<bool U>
inline Integer& Integer::random_lessthan (Integer& r, const unsigned long & m)
{ return Integer::random_lessthan_2exp<U>(r,m);}

template<bool U,class T>
inline Integer Integer::random_lessthan (const T & m)
{
	Integer res ;
	return random_lessthan<U>(res,m);
}

/* ********************************* */
/*  random number of same size as s  */
/* ********************************* */

//! returns a reference to a random number \p r of the size of \p s, exactly.
template<bool U>
inline Integer& Integer::random_exact (Integer& r, const Integer & s)
{
	size_t t = s.bitsize() ;
	random_exact_2exp<U>(r,t);
	return r;
}


/* ************************* */
/*  random number of size m  */
/* ************************* */

//! returns a reference to a random number \p r of the size \p m bits, exactly.
template<bool U>
inline Integer& Integer::random_exact_2exp (Integer& r, const unsigned long int & m)
{
	if (m) random_lessthan_2exp<true>(r,m-1);
	mpz_setbit( (mpz_ptr) &(r.gmp_rep) , m-1);
	if(!U) if (Integer::RandBool()) Integer::negin(r);
	return r;
}

// synonym
template<bool U>
inline Integer& Integer::random_exact (Integer& r, const unsigned long int & m)
{ return Integer::random_exact_2exp<U>(r,m) ; }

template<bool U,class T>
inline Integer Integer::random_exact (const T & s)
{
	Integer res ;
	return random_exact<U>(res,s);
}

/* **************************** */
/*  random number in [[m,M-1]]  */
/* **************************** */

inline Integer& Integer::random_between (Integer& r, const Integer& m, const Integer&M)
{
	assert(M > m);
	random_lessthan(r,Integer(M-m));
	r += m ;
	return (r);
}

inline Integer Integer::random_between (const Integer& m, const Integer &M)
{
	Integer r ;
	return random_between(r,m,M);
}

/* ******************************** */
/*  random number in [[2^m,2^M-1]]  */
/* ******************************** */
// todo : template<bool U, bool V>
inline Integer& Integer::random_between_2exp (Integer& r, const unsigned long int& m, const unsigned long int &M)
{
	assert(M > m);
	r = nonzerorandom((unsigned long int)M-m);
	Integer r1 = random_lessthan_2exp(m);
	r <<= m ;
	r+= r1 ;
	return (r);
}

inline Integer Integer::random_between_2exp (const unsigned long int & m, const unsigned long int &M)
{
	Integer r ;
	return random_between_2exp(r,m,M);
}

// synonym.
inline Integer Integer::random_between (const unsigned long & m, const unsigned long &M)
{ return random_between_2exp(m,M) ; }


inline Integer& Integer::random_between (Integer& r, const unsigned long int& m, const unsigned long int &M)
{
	return random_between_2exp(r,m,M);
}
/* **************/
/*  short hand  */
/* **************/

//! returns a random integer less than...
template<bool U,class T>
inline Integer& Integer::random (Integer& r, const T & m) { return Integer::random_lessthan<U>(r,m) ; }

//! returns a random integer less than...
template<bool U,class T>
inline Integer Integer::random(const T & sz)
{ return Integer::random_lessthan<U,T>(sz); }

inline Integer Integer::random()
{ return Integer::random(sizeof(mp_limb_t)*8) ; }
template<bool U>
inline Integer Integer::random()
{
	Integer rez = Integer::random(sizeof(mp_limb_t)*8) ;
	if (!U) if (Integer::RandBool()) negin(rez);
	return rez;
}

/* *******************/
/*  Non Zero random  */
/* *******************/

template<bool U, class T>
inline Integer Integer::nonzerorandom(const T & sz)
{
	Integer r;
	while(isZero(Integer::random<U,T>(r, sz) )) {} ;
	return r;
}

// BB: It's also 1+random(sz-1)...

template<bool U, class T>
inline Integer& Integer::nonzerorandom (Integer& r, const T& size)
{
	while (isZero(Integer::random<U,T>(r,size))) {} ;
	return r;
}


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
/*  -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*-  */
