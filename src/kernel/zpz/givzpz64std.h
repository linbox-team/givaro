// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givzpz64std.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givzpz64std.h,v 1.15 2011-01-19 18:29:09 bboyer Exp $
// ==========================================================================
// Description:
//   Arithmetic on Z/pZ, with p a prime number less than 2^64
//   Modulo typedef is a signed long number. In case it was modified
//   then bezout algorithm must be changed (coefficient can be negative).
#ifndef _GIVARO_ZPZ64STD_H_
#define _GIVARO_ZPZ64STD_H_
#include "givaro/givinteger.h"
#include "givaro/giverror.h"
#include "givaro/givzpz.h"

// ==========================================================================
// -- This class implement the standard arithmetic with Modulo Elements:
// - The representation of an integer a in Zpz is the value a % p
// ==========================================================================

template<>
class ZpzDom<Std64>
{
	typedef ZpzDom<Std64> Self_t;

public:
	// ----- Exported Types and constantes
	typedef uint64 Residu_t;                    // - type to store residue
	enum { size_rep = sizeof(Residu_t) };      // - size of the storage type
	// ----- Representation of Element of the domain ZpzDom
	typedef int64 Rep;
	typedef int64 Element;


	// ----- Representation of vector of the Element
	typedef Rep* Array;
	typedef const Rep* constArray;

	// ----- Constantes
	const Rep zero;
	const Rep one;

	// ----- Constructor
	ZpzDom() : zero(0), one(1), _p(0) {}
	ZpzDom( Residu_t p, unsigned long e = 1) : zero(0), one(1), _p(p) {}


	Self_t& operator= (const Self_t& D) {
		this->_p = D._p;
		return *this;
	}

	int operator==( const Self_t& BC) const { return _p == BC._p;}
	int operator!=( const Self_t& BC) const { return _p != BC._p;}




	// ----- Access to the modulus
	Residu_t residu() const;
	Residu_t size() const { return _p; }
	Residu_t characteristic() const { return _p; }
	Integer& characteristic(Integer& p) const { return p=_p; }
	Residu_t cardinality() const { return _p; }
	Rep access( const Rep a ) const { return a; }


	// ----- Access to the modulus
	Rep& init( Rep& a ) const;
	void init( const size_t, Array a, constArray b ) const;
	Rep& init( Rep& a, const long i) const ;
	Rep& init( Rep& a, const unsigned long i) const ;
	Rep& init( Rep& a, const long long i) const ;
	Rep& init( Rep& a, const unsigned long long i) const ;
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
	size_t length ( const Rep a ) const;

	// ----- Operations with reduction: r <- a op b mod p, r <- op a mod p
	Rep& mul (Rep& r, const Rep a, const Rep b) const;
	Rep& inv (Rep& r, const Rep a) const;
	Rep& div (Rep& r, const Rep a, const Rep b) const;
	Rep& add (Rep& r, const Rep a, const Rep b) const;
	Rep& sub (Rep& r, const Rep a, const Rep b) const;
	Rep& neg (Rep& r, const Rep a) const;

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
	void axpyin (const size_t sz, Array r, constArray a, constArray x) const;

	// -- axmy: r <- a * x - y mod p
	Rep& axmy  (Rep& r, const Rep a, const Rep b, const Rep c) const;
	void axmy
	(const size_t sz, Array r, constArray a, constArray x, constArray c) const;
	// -- axmyin: r <-  a * x - r  mod p
	Rep& axmyin(Rep& r, const Rep a, const Rep b) const;
	// void axmyin (const size_t sz, Array r, constArray a, constArray x) const;

	// -- maxpy: r <- c - a * b mod p
	Rep& maxpy  (Rep& r, const Rep a, const Rep b, const Rep c) const;
	// -- maxpyin: r <- r - a * b mod p
	Rep& maxpyin(Rep& r, const Rep a, const Rep b) const;
	void maxpyin (const size_t sz, Array r, constArray a, constArray x) const;

	// -- Misc: r <- a mod p
	void assign ( const size_t sz, Array r, constArray a ) const;
	Rep& assign ( Rep& r, const Rep a) const;
	/*
	   Rep& assign ( Rep& r, const long a ) const;
	   Rep& assign ( Rep& r, const unsigned long a ) const;
	   Rep& assign ( Rep& r, const int a ) const;
	   Rep& assign ( Rep& r, const unsigned int a ) const;
	   */
	// ----- random generators
	//     Rep& NONZEROGIVRANDOM(Rep&) const ;
	//     Rep& GIVRANDOM(Rep&) const ;
	template< class RandIter > Rep& random(RandIter&, Rep& r) const ;
	template< class RandIter > Rep& random(RandIter&, Rep& r, long s) const ;
	template< class RandIter > Rep& random(RandIter&, Rep& r, const Rep& b) const ;
	template< class RandIter > Rep& nonzerorandom(RandIter&, Rep& r) const ;
	template< class RandIter > Rep& nonzerorandom(RandIter&, Rep& r, long s) const ;
	template< class RandIter > Rep& nonzerorandom(RandIter&, Rep& r, const Rep& b) const ;


	// <- \sum_i a[i], return 1 if a.size() ==0,
	void reduceadd ( Rep& r, const size_t sz, constArray a ) const;

	// <- \prod_i a[i], return 1 if a.size() ==0,
	void reducemul ( Rep& r, const size_t sz, constArray a ) const;

	// <- \sum_i a[i] * b[i]
	void dotprod ( Rep& r, const size_t sz, constArray a, constArray b ) const;
	void dotprod ( Rep& r, const int bound, const size_t sz, constArray a, constArray b ) const;

	// ----- a -> r: uint64 to double
	void i2d ( const size_t sz, double* r, constArray a ) const;

	// ----- a -> r % p: double to uint64 % p
	void d2i ( const size_t sz, Array r, const double* a ) const;

	// --- IO methods
	std::istream& read ( std::istream& s );
	std::ostream& write( std::ostream& s ) const;
	std::istream& read ( std::istream& s, Rep& a ) const;
	std::ostream& write( std::ostream& s, const Rep a ) const;
	template <class XXX> XXX& convert( XXX& s, const Rep a ) const;
	Integer& write(Integer&, const Rep a ) const;

protected:
	// -- based for modular inverse, d = a*u + b*v
	//   static const int64 gcdext ( int64& u, int64& v, const int64 a, const int64 b );
	int64& gcdext (int64& d, int64& u, int64& v, const int64 a, const int64 b ) const;
	int64& invext (int64& u, const int64 a, const int64 b ) const;

protected:
	// -- data representation of the domain:
	Residu_t _p;

	static void Init();
	static void End();
};


#include "givaro/givzpz64std.inl"
#endif

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
