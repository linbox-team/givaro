// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givzpz16std.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givzpz16std.h,v 1.16 2011-02-02 17:16:43 bboyer Exp $
// ==========================================================================
//
//  Modified by Pascal Giorgi on 2002/02/13  (pascal.giorgi@ens-lyon.fr)
/*! @file givzpz16std.h
 * @ingroup zpz
 * @brief   Arithmetic on Z/pZ, with p a prime number less than 2^14.
 *   Modulo typedef is a signed long number. In case it was modified
 *   then Bézout algorithm must be changed (coefficient can be negative).
 */
#ifndef __GIVARO_zpz16std_H
#define __GIVARO_zpz16std_H


#if 0 /*  Thierry -> JG: Constantes necessaires: */
(typedef) uint16_t : type des int sur 16bits non signe
(typedef) uint32_t : type des int sur 16bits non signe
(typedef) int16_t : type des int sur 16bits signe
(typedef) int32_t : type des int sur 16bits signe
(#define) GIVARO_MAXUINT16: 2^16 -1
(#define) GIVARO_MAXUINT32: 2^32 -1
(#define) GIVARO_MAXULONG: val max en unsigned long
#endif
#if 0
#define GIVARO_BITS_PER_LONGINT        32
#define GIVARO_BITS_PER_INT            32
#define GIVARO_BITS_PER_SHORTINT       16
#define GIVARO_BITS_PER_CHAR           16
typedef char    int8_t;
typedef short   int16_t;
typedef int     int32_t;
typedef unsigned char   uint8_t;
typedef unsigned short  uint16_t;
typedef unsigned int    uint32_t;

#define GIVARO_MAXUINT8                255U            // 2^8-1
#define GIVARO_MAXUINT16               65535U          // 2^16-1
#define GIVARO_MAXUINT32               4294967295U     // 2^32-1
#define GIVARO_MAXULONG                4294967295U     // 2^32-1
#endif
#include "givaro/givbasictype.h"

/*
   Classes d'erreurs:
   GivError::throw_error + des classes d'exception dont les cstos prennent des chaines:
 * GivMathDivZero( " ... " )
 * GivBadFormat( " ... " )
 */
#include "givaro/giverror.h"
// #include "givaro/givzpz16std.h"
#include "givaro/givzpz32std.h"
#include "givaro/giv_randiter.h"

namespace Givaro {

/*! @brief This class implement the standard arithmetic with Modulo Elements.
 * - The representation of an integer a in Zpz is the value a % p
 * - m max is 32768
 * - p max is 32749
 * .
 */
template<>
class ZpzDom<Std16> {
public:
	// ----- Exported Types and constantes
	typedef uint16_t Residu_t;                    // - type to store residue
	enum { size_rep = sizeof(Residu_t) };      // - size of the storage type
	// ----- Representation of Element of the domain ZpzDom
	typedef uint16_t Rep;
	typedef uint16_t Element;
	typedef Element* Element_ptr ;
	typedef const Element* ConstElement_ptr;



	// ----- Representation of vector of the Element
	typedef Rep* Array;
	typedef const Rep* constArray;

	// ----- Constantes
	const Rep zero;
	const Rep one;
	const Rep mOne;

	// ----- Constructor
	ZpzDom() : zero(0), one(1), mOne((Rep)-1), _p(0) {}
	ZpzDom( Residu_t p ) : zero(0), one(1), mOne(Rep(p-1)), _p(p) {}
	ZpzDom( const ZpzDom<Std16>& F) : zero(F.zero), one(F.one), mOne(F.mOne), _p(F._p) {}


	int operator==( const ZpzDom<Std16>& BC) const { return _p == BC._p;}
	int operator!=( const ZpzDom<Std16>& BC) const { return _p != BC._p;}

	ZpzDom<Std16>& operator=( const ZpzDom<Std16>& F)
	{

		F.assign(const_cast<Element&>(one),F.one);
		F.assign(const_cast<Element&>(zero),F.zero);
		F.assign(const_cast<Element&>(mOne),F.mOne);

		this->_p = F._p;
		return *this;
	}

	// ----- Access to the modulus
	Residu_t residu() const;
	Residu_t size() const { return _p; }
	Residu_t characteristic() const { return _p; }
	Integer& characteristic( Integer& p) const { return p=_p; }
	Rep access( const Rep a ) const { return a; }

	// ----- Convert from Element to int
	int16_t& convert( int16_t& x , const Rep a) const { return x=(int16_t)(a);}
	uint16_t& convert( uint16_t& x , const Rep a) const { return x=(uint16_t)(a);}
	unsigned long& convert( unsigned long& x , const Rep a) const { return x=(unsigned long)(a);}
	double& convert( double& x , const Rep a) const { return x=(double)(a);}
	int& convert( int& x , const Rep a) const { return x=int(a);}
	Integer& convert(Integer& i, const Rep a) const {
		unsigned long ur;
		return i = (Integer)convert(ur, a);
	}



	// ----- Access to the modulus
	Rep& init( Rep& a ) const;
	void init( const size_t, Array a, constArray b ) const;
	Rep& init( Rep& a, const long i) const ;
	Rep& init( Rep& a, const unsigned long i) const ;
	Rep& init( Rep& a, const int i) const ;
	Rep& init( Rep& a, const unsigned int i) const ;
	Rep& init( Rep& a, const double i) const ;
	Rep& init( Rep& a, const float i) const ;
	Rep& init( Rep& a, const Integer& i) const ;

	// ----- Misc methods
	int areEqual( const  Rep, const Rep) const;
	int areNEqual( const Rep, const Rep) const;
	int isZero( const Rep a ) const;
	int isnzero( const Rep a ) const;
	int isOne ( const Rep a ) const;
	int isMOne ( const Rep a ) const;
	size_t length ( const Rep a ) const;

	// ----- Operations with reduction: r <- a op b mod p, r <- op a mod p
	Rep& mul (Rep& r, const Rep a, const Rep b) const;
	Rep& div (Rep& r, const Rep a, const Rep b) const;
	Rep& add (Rep& r, const Rep a, const Rep b) const;
	Rep& sub (Rep& r, const Rep a, const Rep b) const;
	Rep& neg (Rep& r, const Rep a) const;
	Rep& inv (Rep& r, const Rep a) const;

	Rep& mulin (Rep& r, const Rep a) const;
	Rep& divin (Rep& r, const Rep a) const;
	Rep& addin (Rep& r, const Rep a) const;
	Rep& subin (Rep& r, const Rep a) const;
	Rep& negin (Rep& r) const;
	Rep& invin (Rep& r) const;

	// ----- Operations with reduction: r <- a op b mod p, r <- op a mod p
	void mul (const size_t sz, Array r, constArray a, constArray b) const;
	void mul (const size_t sz, Array r, constArray a, Rep b) const;

	void div (const size_t sz, Array r, constArray a, constArray b) const;
	void div (const size_t sz, Array r, constArray a, Rep b) const;

	void add (const size_t sz, Array r, constArray a, constArray b) const;
	void add (const size_t sz, Array r, constArray a, Rep b) const;

	void sub (const size_t sz, Array r, constArray a, constArray b) const;
	void sub (const size_t sz, Array r, constArray a, Rep b) const;

	void neg (const size_t sz, Array r, constArray a) const;
	void inv (const size_t sz, Array r, constArray a) const;

	// -- axpy: r <- a * x + y mod p
	Rep& axpy  (Rep& r, const Rep a, const Rep b, const Rep c) const;
	void axpy
	(const size_t sz, Array r, constArray a, constArray x, constArray c) const;
	// -- axpyin: r <- r + a * x mod p
	Rep& axpyin(Rep& r, const Rep a, const Rep b) const;
	void axpyin
	(const size_t sz, Array r, constArray a, constArray x) const;

	// -- maxpy: r <- c - a * b mod p
	Rep& maxpy (Rep& r, const Rep a, const Rep b, const Rep c) const;
	// -- maxpyin: r <- r - a * x mod p
	Rep& maxpyin (Rep& r, const Rep a, const Rep b) const;
	void maxpyin (const size_t sz, Array r, constArray a, constArray x) const;

	// -- axmy: r <- a * x - y mod p
	Rep& axmy  (Rep& r, const Rep a, const Rep b, const Rep c) const;
	void axmy
	(const size_t sz, Array r, constArray a, constArray x, constArray c) const;
	// -- axmyin: r <-  a * x - r mod p
	Rep& axmyin(Rep& r, const Rep a, const Rep b) const;
	// void axmyin (const size_t sz, Array r, constArray a, constArray x) const;

	// -- Misc: r <- a mod p
	void assign ( const size_t sz, Array r, constArray a ) const;
    	Rep& assign ( Rep& r, const Rep a) const;

	// ----- random generators
	template< class RandIter > Rep& random(RandIter&, Rep& r) const ;
	template< class RandIter > Rep& random(RandIter&, Rep& r, long s) const ;
	template< class RandIter > Rep& random(RandIter&, Rep& r, const Rep& b) const ;
	template< class RandIter > Rep& nonzerorandom(RandIter&, Rep& r) const ;
	template< class RandIter > Rep& nonzerorandom(RandIter&, Rep& r, long s) const ;
	template< class RandIter > Rep& nonzerorandom(RandIter&, Rep& r, const Rep& b) const ;

	typedef GIV_randIter< ZpzDom<Std16>, Rep > randIter;

	// <- \sum_i a[i], return 1 if a.size() ==0,
	void reduceadd ( Rep& r, const size_t sz, constArray a ) const;

	// <- \prod_i a[i], return 1 if a.size() ==0,
	void reducemul ( Rep& r, const size_t sz, constArray a ) const;

	// <- \sum_i a[i] * b[i]
	void dotprod ( Rep& r, const size_t sz, constArray a, constArray b ) const;
	void dotprod ( Rep& r, const int bound, const size_t sz, constArray a, constArray b ) const;

	// ----- a -> r: uint16_t to double
	void i2d ( const size_t sz, double* r, constArray a ) const;

	// ----- a -> r % p: double to uint16_t % p
	void d2i ( const size_t sz, Array r, const double* a ) const;

	// --- IO methods
	std::istream& read ( std::istream& s );
	std::ostream& write( std::ostream& s ) const;
	std::istream& read ( std::istream& s, Rep& a ) const;
	std::ostream& write( std::ostream& s, const Rep a ) const;

protected:
	// -- based for modular inverse, d = a*u + b*v
	//   static const int32_t gcdext ( int32_t& u, int32_t& v, const int32_t a, const int32_t b );
	int32_t& gcdext (int32_t& d, int32_t& u, int32_t& v, const int32_t a, const int32_t b ) const;
	int32_t& invext (int32_t& u, const int32_t a, const int32_t b ) const;

protected:
	// -- data representation of the domain:
	Residu_t _p;

	static void Init();
	static void End();

public: static inline Residu_t getMaxModulus() { return 32768; }

};

} // namespace Givaro

#include "givaro/givzpz16std.inl"

#endif // __GIVARO_zpz16std_H
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
