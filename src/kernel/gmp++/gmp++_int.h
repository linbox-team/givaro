// ========================================================================
// Copyright(c)'1994-2010 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: M. Samama, T. Gautier, JG. Dumas
// Modified by: B. Boyer
// Time-stamp: <22 Oct 10 15:35:39 Jean-Guillaume.Dumas@imag.fr>
// ========================================================================
// Description:
// Integer class definition based on Gmp (>V2.0 or 1.3.2)
#ifndef __GIVARO_GMPplusplus_integer_H
#define __GIVARO_GMPplusplus_integer_H

/*! @file kernel/gmp++/gmp++_int.h
 * @ingroup integers
 * Core gmp++_int.h.
 */

#include <vector>
#include <list>
#include <string>
#include <assert.h>
// #include <iostream>

#ifndef __GIVARO_GMPplusplus_H
#warning "you should include <gmp++/gmp++.h> before <gmp++/gmp++_int.h> (or prepare for the worse)"
#endif

#ifdef __USE_64_bits__
#    define __USE_GMPPLUSPLUS_SIXTYFOUR__
#endif

#ifdef __USE_ISOC99
#    define __USE_GMPPLUSPLUS_SIXTYFOUR__
#endif

namespace Givaro {

	//------------------------------------------------------ Friend Integer
	// forward declaration.
	class Integer;

	// (FILE gmp++_int_gcd.C)
	Integer& 	inv (Integer& u, const Integer& a, const Integer& b);
	Integer& 	invin (Integer& u, const Integer& b);

	Integer 	gcd (const Integer& a, const Integer& b);
	Integer 	gcd (Integer& u, Integer& v,const Integer& a, const Integer& b);
	Integer& 	gcd (Integer& g, const Integer& a, const Integer& b);
	Integer& 	gcd (Integer& g, Integer& u, Integer& v, const Integer& a, const Integer& b);
	Integer 	pp( const Integer& P, const Integer& Q );
	Integer& 	lcm (Integer& g, const Integer& a, const Integer& b);
	Integer 	lcm (const Integer& a, const Integer& b);

	// (FILE gmp++_int_pow.C)
	Integer& 	pow(Integer& Res, const Integer& n, const long int l);
	Integer& 	pow(Integer& Res, const long unsigned int n, const long unsigned int l);
	Integer& 	pow(Integer& Res, const Integer& n, const long unsigned int l);
	Integer& 	pow(Integer& Res, const Integer& n, const int l) ;
	Integer& 	pow(Integer& Res, const Integer& n, const unsigned int l) ;
	Integer 	pow(const Integer& n, const long int l);
	Integer 	pow(const Integer& n, const long unsigned int l);
	Integer 	pow(const Integer& n, const int l) ;
	Integer 	pow(const Integer& n, const unsigned int l);

	Integer& 	powmod(Integer& Res, const Integer& n, const long unsigned int e, const Integer& m);
	Integer& 	powmod(Integer& Res, const Integer& n, const long int e, const Integer& m);
	Integer& 	powmod(Integer& Res, const Integer& n, const unsigned int e, const Integer& m) ;
	Integer& 	powmod(Integer& Res, const Integer& n, const int e, const Integer& m)  ;
	Integer& 	powmod(Integer& Res, const Integer& n, const Integer& e, const Integer& m);
	Integer 	powmod(const Integer& n, const long unsigned int e, const Integer& m);
	Integer 	powmod(const Integer& n, const long int e, const Integer& m);
	Integer 	powmod(const Integer& n, const unsigned int e, const Integer& m) ;
	Integer 	powmod(const Integer& n, const int e, const Integer& m) ;
	Integer 	powmod(const Integer& n, const Integer& e, const Integer& m);

	// (FILE inline)
	int 		sign   (const Integer& a);

	// (FILE gmp++_int_compare.C)
	int 		compare(const Integer& a, const Integer& b);

	int 		absCompare(const Integer& a, const Integer& b);
	int 		absCompare(const Integer& a, const double b);
	int 		absCompare(const Integer& a, const float b);
	int 		absCompare(const Integer& a, const long unsigned b);
	int 		absCompare(const Integer& a, const unsigned b);
	int 		absCompare(const Integer& a, const long int b);
	int 		absCompare(const Integer& a, const int b);
	template<class T>
	int             absCompare( const T a,  const Integer & b)
	{
		return absCompare(b,a);
	}

	int 		isZero (const Integer& a);
	int 		nonZero (const Integer& a);
	int 		isOne  (const Integer& a);

	// (FILE gmp++_int_misc.C)
	Integer 	fact ( long unsigned int l);

	Integer 	sqrt(const Integer& p);
	Integer 	sqrtrem(const Integer& p, Integer& rem);
	Integer& 	sqrt(Integer& r, const Integer& p);
	Integer& 	sqrtrem(Integer& r, const Integer& p, Integer& rem);
	bool 		root(Integer& q, const Integer&, unsigned int n);

	long 		logp(const Integer& a, const Integer& p) ;
	double 		logtwo(const Integer& a) ;
	double 		naturallog(const Integer& a) ;

	void 		swap(Integer& , Integer&);

	int 		isperfectpower  (const Integer& );

	Integer 	abs(const Integer& n);

	Integer& 	prevprime(Integer&, const Integer& p);
	Integer& 	nextprime(Integer&, const Integer& p);
	int 		probab_prime(const Integer& p);
	int 		probab_prime(const Integer& p, int r);
	int 		jacobi(const Integer& u, const Integer& v) ;
	int 		legendre(const Integer& u, const Integer& v) ;

	long unsigned 	length (const Integer& a);
	// (FILE gmp++_int_io.C)

	std::istream& 	operator >> (std::istream &i, Integer& n);
	std::ostream& 	operator << (std::ostream &o, const Integer& n);
	std::ostream& 	absOutput (std::ostream &o, const Integer& n);
	void 		importWords(Integer&, size_t, int, int, int, size_t, const void*);

	//------------------------------------------------------ Class Integer
	/*! @ingroup integers
	 * This is the Integer class.
	 * An Integer is represented as a GMP integer.
	 * This class provides arithmetic on Integers.
	 */
	class Integer {

	public:
		//! vector of limbs (ie a gmp number).
		typedef std::vector<mp_limb_t> vect_t;
		//--------------------------------------cstors & dstors
		// (FILE gmp++_cstor.C)
		/*! @name Constructor/Destructors.
		 * Constructors and destructor for an Integer.
		*/
		/// Constructor form a known type.
		///@param n input to be constructed from
		giv_all_inlined Integer(int n = 0);
		//! @overload Givaro::Integer(int)
		giv_all_inlined Integer(long int n);
		//! @overload Givaro::Integer(int)
		giv_all_inlined Integer(unsigned char n);
		//! @overload Givaro::Integer(int)
		giv_all_inlined Integer(unsigned int n);
		//! @overload Givaro::Integer(int)
		giv_all_inlined Integer(long unsigned int n);
#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
		//! @overload Givaro::Integer(int)
		giv_all_inlined Integer(long long int n);
		//! @overload Givaro::Integer(int)
		giv_all_inlined Integer(long long unsigned int n);
#endif
		//! @overload Givaro::Integer(int)
		giv_all_inlined Integer(double n);
		//! @overload Givaro::Integer(int)
		giv_all_inlined Integer(const char *n);

		/*! Copy constructor
		 * @param n input to be constructed from
		 */
		giv_all_inlined Integer(const Integer& n);

		/*! Creates a new Integer from pointers.
		 * @param d array
		 * @param sz size
		 */
		giv_all_inlined Integer(long unsigned* d, long sz);
		/*! Creates a new Integers for a vector of limbs
		* @param v vector of limbs
		*/
		giv_all_inlined Integer( const vect_t &v);

		//! destructor
		giv_all_inlined ~Integer();
		///@}

		//------------------------------------- predefined null and one
		//! zero (0)
		static const Integer zero ;
		//! one (1)
		static const Integer one ;
		//! minus one (-1)
		static const Integer mOne ;



		// -- Assignment and copy operators
		/*! @name Assignment and copy operators */
		/*! copy from an integer.
		 * @param n integer to copy.
		 */
		///@{
		giv_all_inlined Integer& operator = (const Integer& n);
		giv_all_inlined Integer& logcpy(const Integer& n);
		giv_all_inlined Integer& copy(const Integer& n);
		///@}

		//------------------Equalities and inequalities between integers and longs
		// (FILE gmp++_int_compare.C)
		//! @name Comparisons functions.
		///@{
		/*! Compares two integers.
		 * @param a integer
		 * @param b integer
		 * @return \c 1 if \f$a > b\f$, \c 0 if \f$a = b\f$ and \p -1 otherwise.
		 */
		giv_all_inlined friend int compare(const Integer& a, const Integer& b);

		/** Compare the norm of two integers.
		 * @param a integer
		 * @param b integer
		 * @return \c 1 if \f$|a| > |b|\f$, \c 0 if \f$|a| = |b|\f$ and \p -1 otherwise.
		 */
		giv_all_inlined friend int absCompare(const Integer& a, const Integer& b);
		/** @overload Integer::absCompare(Integer, Integer) */
		giv_all_inlined friend int absCompare(const Integer& a, const double b);
		/** @overload Integer::absCompare(Integer, Integer) */
		giv_all_inlined friend int absCompare(const Integer& a, const float b);
		/** @overload Integer::absCompare(Integer, Integer) */
		giv_all_inlined friend int absCompare(const Integer& a, const long unsigned b);
		/** @overload Integer::absCompare(Integer, Integer) */
		giv_all_inlined friend int absCompare(const Integer& a, const unsigned b);
		/** @overload Integer::absCompare(Integer, Integer) */
		giv_all_inlined friend int absCompare(const Integer& a, const long int b);
		/** @overload Integer::absCompare(Integer, Integer) */
		giv_all_inlined friend int absCompare(const Integer& a, const int b);
		/** @overload Integer::absCompare(Integer, Integer) */
		template<class T>
		giv_all_inlined friend int absCompare( const T a,  const Integer & b) ;

		//! name compare to 1 and 0
		//! @param a
		friend giv_all_inlined  int isOne(const Integer& a);



			//! name compare to 1 and 0
		//! @param a
friend giv_all_inlined  int nonZero(const Integer& a);

		//! name compare to 1 and 0
		//! @param a
	friend giv_all_inlined  int isZero(const Integer& a);
	//! @overload Givaro::isZero(Integer);
		friend giv_all_inlined  int isZero(const short int a);
	//! @overload Givaro::isZero(Integer);
		friend giv_all_inlined  int isZero(const int a);
	//! @overload Givaro::isZero(Integer);
		friend giv_all_inlined  int isZero(const long int a);
	//! @overload Givaro::isZero(Integer);
		friend giv_all_inlined  int isZero(const unsigned short int a);
	//! @overload Givaro::isZero(Integer);
		friend giv_all_inlined  int isZero(const unsigned int a);
	//! @overload Givaro::isZero(Integer);
		friend giv_all_inlined  int isZero(const long unsigned int a);
#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
#if 1 /*  use of C++0x long long integer constant */
	//! @overload Givaro::isZero(Integer);
		friend giv_all_inlined  int isZero(const long long unsigned int a);
#endif
	//! @overload Givaro::isZero(Integer);
		friend giv_all_inlined  int isZero(const long long int a);
#endif
		//! isleq
		//! @param a,b
		template<class A,class B>
		static giv_all_inlined bool isleq(const A&a,const B&b)
		{
		       	return a<=b ;
	       	}

		///@}


		//! @name Comparison operators.
		///@{
		//! greater or equal
		/// @param l integer to be compared to
		giv_all_inlined int operator >= (const Integer & l) const;
		/** @overload Integer::operator>=(Integer) */
		giv_all_inlined int operator >= (const int l) const;
		/** @overload Integer::operator>=(Integer) */
		giv_all_inlined int operator >= (const long int l) const;
		/** @overload Integer::operator>=(Integer) */
		giv_all_inlined int operator >= (const long unsigned int l) const;
		/** @overload Integer::operator>=(Integer) */
		giv_all_inlined int operator >= (const unsigned int l) const;
		/** @overload Integer::operator>=(Integer) */
		giv_all_inlined int operator >= (const double l) const;
		/** @overload Integer::operator>=(Integer) */
		giv_all_inlined int operator >= (const float l) const;

		//! greater or equal.
		/// @param l,n integers to compare
		giv_all_inlined friend int operator >= (unsigned int l, const Integer& n);
		/** @overload Integer::operator>=(unsigned, Integer) */
		giv_all_inlined friend int operator >= (float l, const Integer& n);
		/** @overload Integer::operator>=(unsigned, Integer) */
		giv_all_inlined friend int operator >= (double l, const Integer& n);
		/** @overload Integer::operator>=(unsigned, Integer) */
		giv_all_inlined friend int operator >= (int l, const Integer& n);
		/** @overload Integer::operator>=(unsigned, Integer) */
		giv_all_inlined friend int operator >= (long int l, const Integer& n);
		/** @overload Integer::operator>=(unsigned, Integer) */
		giv_all_inlined friend int operator >= (long unsigned int l, const Integer& n);


		//! less or equal
		/// @param l integer to be compared to
		giv_all_inlined int operator <= (const Integer & l) const;
		/** @overload Integer::operator<=(Integer) */
		giv_all_inlined int operator <= (const int l) const;
		/** @overload Integer::operator<=(Integer) */
		giv_all_inlined int operator <= (const long int l) const;
		/** @overload Integer::operator<=(Integer) */
		giv_all_inlined int operator <= (const long unsigned int l) const;
		/** @overload Integer::operator<=(Integer) */
		giv_all_inlined int operator <= (const unsigned int l) const;
		/** @overload Integer::operator<=(Integer) */
		giv_all_inlined int operator <= (const double l) const;
		/** @overload Integer::operator<=(Integer) */
		giv_all_inlined int operator <= (const float l) const;

		//! less or equal
		/// @param l,n integers to compare
		giv_all_inlined friend int operator <= (unsigned int l, const Integer& n);
		/** @overload Integer::operator>=(unsigned, Integer) */
		giv_all_inlined friend int operator <= (float l, const Integer& n);
		/** @overload Integer::operator>=(unsigned, Integer) */
		giv_all_inlined friend int operator <= (double l, const Integer& n);
		/** @overload Integer::operator>=(unsigned, Integer) */
		giv_all_inlined friend int operator <= (int l, const Integer& n);
		/** @overload Integer::operator>=(unsigned, Integer) */
		giv_all_inlined friend int operator <= (long int l, const Integer& n);
		/** @overload Integer::operator>=(unsigned, Integer) */
		giv_all_inlined friend int operator <= (long unsigned int l, const Integer& n);

		/*! operator != (not equal)
		 * @param l integer
		 * @return \c 1 iff l == this
		 */
		giv_all_inlined int operator != (const Integer & l) const;
		/** @overload Integer::operator!=(Integer) */
		giv_all_inlined int operator != (const int l) const;
		/** @overload Integer::operator!=(Integer) */
		giv_all_inlined int operator != (const long int l) const;
		/** @overload Integer::operator!=(Integer) */
		giv_all_inlined int operator != (const long unsigned int l) const;
		/** @overload Integer::operator!=(Integer) */
		giv_all_inlined int operator != (const unsigned int l) const;
		/** @overload Integer::operator!=(Integer) */
		giv_all_inlined int operator != (const double l) const;
		/** @overload Integer::operator!=(Integer) */
		giv_all_inlined int operator != (const float l) const;

		/*! operator != (not equal)
		 * @param l,n integer
		 * @return \c 1 iff l == n
		 */
		giv_all_inlined friend int operator != (unsigned int l, const Integer& n);
		/** @overload Integer::operator!=(unsigned, Integer) */
		giv_all_inlined friend int operator != (float l, const Integer& n);
		/** @overload Integer::operator!=(unsigned, Integer) */
		giv_all_inlined friend int operator != (double l, const Integer& n);
		/** @overload Integer::operator!=(unsigned, Integer) */
		giv_all_inlined friend int operator != (int l, const Integer& n);
		/** @overload Integer::operator!=(unsigned, Integer) */
		giv_all_inlined friend int operator != (long int l, const Integer& n);
		/** @overload Integer::operator!=(unsigned, Integer) */
		giv_all_inlined friend int operator != (long unsigned int l, const Integer& n);


		//! Equality
		/// @param l integer to be compared to
		giv_all_inlined int operator == (const Integer & l) const;
		/** @overload Integer::operator==(Integer) */
		giv_all_inlined int operator == (const int l) const;
		/** @overload Integer::operator==(Integer) */
		giv_all_inlined int operator == (const long int l) const;
		/** @overload Integer::operator==(Integer) */
		giv_all_inlined int operator == (const long unsigned int l) const;
		/** @overload Integer::operator==(Integer) */
		giv_all_inlined int operator == (const unsigned int l) const;
		/** @overload Integer::operator==(Integer) */
		giv_all_inlined int operator == (const double l) const;
		/** @overload Integer::operator==(Integer) */
		giv_all_inlined int operator == (const float l) const;

		//! Equality
		/// @param l,n integers to compare
		giv_all_inlined friend int operator == (unsigned int l, const Integer& n);
		/** @overload Integer::operator==(unsigned, Integer) */
		giv_all_inlined friend int operator == (float l, const Integer& n);
		/** @overload Integer::operator==(unsigned, Integer) */
		giv_all_inlined friend int operator == (double l, const Integer& n);
		/** @overload Integer::operator==(unsigned, Integer) */
		giv_all_inlined friend int operator == (int l, const Integer& n);
		/** @overload Integer::operator==(unsigned, Integer) */
		giv_all_inlined friend int operator == (long int l, const Integer& n);
		/** @overload Integer::operator==(unsigned, Integer) */
		giv_all_inlined friend int operator == (long unsigned int l, const Integer& n);


		//! greater (strict)
		/// @param l integer to be compared to
		giv_all_inlined int operator > (const Integer & l) const;
		/** @overload Integer::operator>(Integer) */
		giv_all_inlined int operator > (const int l) const;
		/** @overload Integer::operator>(Integer) */
		giv_all_inlined int operator > (const long int l) const;
		/** @overload Integer::operator>(Integer) */
		giv_all_inlined int operator > (const long unsigned int l) const;
		/** @overload Integer::operator>(Integer) */
		giv_all_inlined int operator > (const unsigned int l) const;
		/** @overload Integer::operator>(Integer) */
		giv_all_inlined int operator > (const double l) const;
		/** @overload Integer::operator>(Integer) */
		giv_all_inlined int operator > (const float l) const;

		//! greater (strict)
		/// @param l,n integers to compare
		giv_all_inlined friend int operator > (unsigned int l, const Integer& n);
		/** @overload Integer::operator>(unsigned, Integer) */
		giv_all_inlined friend int operator > (float l, const Integer& n);
		/** @overload Integer::operator>(unsigned, Integer) */
		giv_all_inlined friend int operator > (double l, const Integer& n);
		/** @overload Integer::operator>(unsigned, Integer) */
		giv_all_inlined friend int operator > (int l, const Integer& n);
		/** @overload Integer::operator>(unsigned, Integer) */
		giv_all_inlined friend int operator > (long int l, const Integer& n);
		/** @overload Integer::operator>(unsigned, Integer) */
		giv_all_inlined friend int operator > (long unsigned int l, const Integer& n);

		//! less (strict)
		/// @param l integer to be compared to
		giv_all_inlined int operator < (const Integer & l) const;
		/** @overload Integer::operator<(Integer) */
		giv_all_inlined int operator < (const int l) const;
		/** @overload Integer::operator<(Integer) */
		giv_all_inlined int operator < (const long int l) const;
		/** @overload Integer::operator<(Integer) */
		giv_all_inlined int operator < (const long unsigned int l) const;
		/** @overload Integer::operator<(Integer) */
		giv_all_inlined int operator < (const unsigned int l) const;
		/** @overload Integer::operator<(Integer) */
		giv_all_inlined int operator < (const double l) const;
		/** @overload Integer::operator<(Integer) */
		giv_all_inlined int operator < (const float l) const;

		//! less (strict)
		/// @param l,n integers to compare
		giv_all_inlined friend int operator < (unsigned int l, const Integer& n);
		/** @overload Integer::operator<(unsigned, Integer) */
		giv_all_inlined friend int operator < (float l, const Integer& n);
		/** @overload Integer::operator<(unsigned, Integer) */
		giv_all_inlined friend int operator < (double l, const Integer& n);
		/** @overload Integer::operator<(unsigned, Integer) */
		giv_all_inlined friend int operator < (int l, const Integer& n);
		/** @overload Integer::operator<(unsigned, Integer) */
		giv_all_inlined friend int operator < (long int l, const Integer& n);
		/** @overload Integer::operator<(unsigned, Integer) */
		giv_all_inlined friend int operator < (long unsigned int l, const Integer& n);
		///@}

		// ---------------- Bit logic
		// (FILE gmp++_int_misc.C)

		/*! @name Bit logic operators.  */
		///@{
		//! @brief XOR (^)
		//! @param a integer
		giv_all_inlined Integer operator^ (const Integer& a) const;   // XOR
		/** @overload Integer::operator^(Integer) */
		giv_all_inlined Integer operator^ (const long unsigned int & a) const;
		/** @overload Integer::operator^(Integer) */
		giv_all_inlined Integer operator^ (const unsigned int& a) const;
		//! @brief XOR inplace (^=)
		//! @param a integer
		giv_all_inlined Integer& operator^= (const Integer&a);   // XOR
		/** @overload Integer::operator^=(Integer) */
		giv_all_inlined Integer& operator^= (const long unsigned int & a);   // XOR
		/** @overload Integer::operator^=(Integer) */
		giv_all_inlined Integer& operator^= (const unsigned int&a);   // XOR

		//! OR (|)
		//! @param a integer
		giv_all_inlined Integer operator| (const Integer& a ) const;   // OR
		/** @overload Integer::operator|(Integer) */
		giv_all_inlined Integer operator| (const long unsigned int & a) const;
		/** @overload Integer::operator|(Integer) */
		giv_all_inlined Integer operator| (const unsigned int& a) const;
		//! OR inplace (|=)
		//! @param a integer
		giv_all_inlined Integer& operator|= (const Integer& a );   // OR
		/** @overload Integer::operator|=(Integer) */
		giv_all_inlined Integer& operator|= (const long unsigned int & a );   // OR
		/** @overload Integer::operator|=(Integer) */
		giv_all_inlined Integer& operator|= (const unsigned int& a );   // OR

		//! AND (&)
		//! @param a integer
		giv_all_inlined Integer operator& (const Integer&a) const;   // AND
		/** @overload Integer::operator&(Integer) */
		giv_all_inlined unsigned int operator& (const unsigned int& a) const;
		/** @overload Integer::operator&(Integer) */
		giv_all_inlined long unsigned operator& (const long unsigned int & a) const;

		//! AND inplace (&=)
		//! @param a integer
		giv_all_inlined Integer& operator&= (const Integer&a);   // AND
		/** @overload Integer::operator&=(Integer) */
		giv_all_inlined Integer& operator&= (const long unsigned int &a);   // AND
		/** @overload Integer::operator&=(Integer) */
		giv_all_inlined Integer& operator&= (const unsigned int&a);   // AND

		//! complement to 1 (~)
		giv_all_inlined Integer operator ~ () const;   // 1 complement

		//! left shift (<<)
		//! @param l shift
		giv_all_inlined Integer operator<< (int l) const; // lshift
		/** @overload Integer::operator<<(int) */
		giv_all_inlined Integer operator<< (long int l) const; // lshift
		/** @overload Integer::operator<<(int) */
		giv_all_inlined Integer operator<< (unsigned int l) const; // lshift
		/** @overload Integer::operator<<(int) */
		giv_all_inlined Integer operator<< (long unsigned int l) const; // lshift

		//! left shift inplace (<<=)
		//! @param l shift
		giv_all_inlined Integer& operator<<= (int l) ; // lshift
		/** @overload Integer::operator<<=(int) */
		giv_all_inlined Integer& operator<<= (long int l) ; // lshift
		/** @overload Integer::operator<<=(Integer) */
		giv_all_inlined Integer& operator<<= (unsigned int l) ; // lshift
		/** @overload Integer::operator<<=(Integer) */
		giv_all_inlined Integer& operator<<= (long unsigned int l) ; // lshift

		//! right shift (>>)
		//! @param l shift
		giv_all_inlined Integer operator>> (int l) const; // rshift
		/** @overload Integer::operator>>(int) */
		giv_all_inlined Integer operator>> (long int l) const; // rshift
		/** @overload Integer::operator>>(int) */
		giv_all_inlined Integer operator>> (unsigned int l) const; // rshift
		/** @overload Integer::operator>>(int) */
		giv_all_inlined Integer operator>> (long unsigned int l) const; // rshift

		//! right shift inplace (>>=)
		//! @param l shift
		giv_all_inlined Integer& operator>>= (int l) ; // rshift
		/** @overload Integer::operator>>=(int) */
		giv_all_inlined Integer& operator>>= (long int l) ; // rshift
		/** @overload Integer::operator>>=(int) */
		giv_all_inlined Integer& operator>>= (unsigned int l) ; // rshift
		/** @overload Integer::operator>>=(int) */
		giv_all_inlined Integer& operator>>= (long unsigned int l) ; // rshift
		///@}


		// - Methods
		/*! @name Addition, substraction, multiplication */
		///@{
		// (FILE gmp++_int_add.C)
		/*!  Addition (inplace)
		 * <code>res+=n</code>.
		 * @param res as in the formula
		 * @param n as in the formula
		 */
		static giv_all_inlined  Integer& addin (Integer& res, const Integer& n);
		/** @overload Integer::addin(Integer,Integer) */
		static giv_all_inlined  Integer& addin (Integer& res, const long int n);
		/** @overload Integer::addin(Integer,Integer) */
		static giv_all_inlined  Integer& addin (Integer& res, const long unsigned int n);

		/*!  Addition
		 * <code>res=n1+n2</code>.
		 * @param res as in the formula
		 * @param n1 as in the formula
		 * @param n2 as in the formula
		 */
		static giv_all_inlined  Integer& add   (Integer& res, const Integer& n1, const Integer& n2);
		/** @overload Integer::add(Integer,Integer,Integer) */
		static giv_all_inlined  Integer& add   (Integer& res, const Integer& n1, const long int n2);
		/** @overload Integer::add(Integer,Integer,Integer) */
		static giv_all_inlined  Integer& add   (Integer& res, const Integer& n1, const long unsigned int n2);

		// (FILE gmp++_int_sub.C)
		/*!  Substraction (inplace)
		 * <code>res-=n</code>.
		 * @param res as in the formula
		 * @param n as in the formula
		 */
		static giv_all_inlined  Integer& subin (Integer& res, const Integer& n);
		/** @overload Integer::subin(Integer,Integer) */
		static giv_all_inlined  Integer& subin (Integer& res, const long int n);
		/** @overload Integer::subin(Integer,Integer) */
		static giv_all_inlined  Integer& subin (Integer& res, const long unsigned int n);

		/*!  Substraction
		 * <code>res=n1-n2</code>.
		 * @param res as in the formula
		 * @param n1 as in the formula
		 * @param n2 as in the formula
		 */
		static giv_all_inlined  Integer& sub   (Integer& res, const Integer& n1, const Integer& n2);
		/** @overload Integer::sub(Integer,Integer,Integer) */
		static giv_all_inlined  Integer& sub   (Integer& res, const Integer& n1, const long int n2);
		/** @overload Integer::sub(Integer,Integer,Integer) */
		static giv_all_inlined  Integer& sub   (Integer& res, const Integer& n1, const long unsigned int n2);

		/*!  Negation (inplace)
		 * <code>res=-res</code>.
		 * @param res as in the formula
		 */
		static giv_all_inlined  Integer& negin (Integer& res);
		/*!  Negation
		 * <code>res=-n</code>.
		 * @param n as in the formula
		 * @param res as in the formula
		 */
		static giv_all_inlined  Integer& neg   (Integer& res, const Integer& n);

		/*!  Multiplication (inplace)
		 * <code>res*=n</code>.
		 * @param res as in the formula
		 * @param n as in the formula
		 */
		static giv_all_inlined  Integer& mulin (Integer& res, const Integer& n);
		/** @overload Integer::mulin(Integer,Integer) */
		static giv_all_inlined  Integer& mulin (Integer& res, const long int n);
		/** @overload Integer::mulin(Integer,Integer) */
		static giv_all_inlined  Integer& mulin (Integer& res, const long unsigned int n);

		/*! Multiplication
		 * <code>res=n1*n2</code>.
		 * @param res as in the formula
		 * @param n1 as in the formula
		 * @param n2 as in the formula
		 */
		static giv_all_inlined  Integer& mul   (Integer& res, const Integer& n1, const Integer& n2);
		/** @overload Integer::mul(Integer,Integer,Integer) */
		static giv_all_inlined  Integer& mul   (Integer& res, const Integer& n1, const long int n2);
		/** @overload Integer::mul(Integer,Integer,Integer) */
		static giv_all_inlined  Integer& mul   (Integer& res, const Integer& n1, const long unsigned int n2);
		///@}

		//----------------Elementary arithmetic between Integers & longs
		// (FILE gmp++_int_add.C)
		/*! @name Addition, substraction, multiplication operators*/
		///@{
		/*! operator \c +.
		 * @return <code> (*this)+n</code>
		 * @param n as in the formula.
		 */
		giv_all_inlined Integer  operator + (const Integer& n) const;
		/** @overload Integer::operator+(Integer) */
		giv_all_inlined Integer  operator + (const long unsigned int n) const;
		/** @overload Integer::operator+(Integer) */
		giv_all_inlined Integer  operator + (const long int n) const;

		/*! operator \c += .
		 * @param n asfriend In the formula.
		 * @return <code> (*this) += n</code>.
		 */
		giv_all_inlined Integer& operator += (const Integer& n);
		/** @overload Integer::operator+=(Integer) */
		giv_all_inlined Integer& operator += (const long unsigned int n);
		/** @overload Integer::operator+=(Integer) */
		giv_all_inlined Integer& operator += (const long int n);
		/** @overload Integer::operator+=(Integer) */
		template<class XXX>
		Integer& operator +=(const XXX& n) {
			return this->operator += ( (Integer)n );
		}

		// - Friends
		//! operator +.
		//! @param l,n to be added
		friend giv_all_inlined  Integer operator + (const int l, const Integer& n);
		/** @overload friend Integer::operator+(int,Integer) */
		friend giv_all_inlined  Integer operator + (const unsigned int l, const Integer& n);
		/** @overload friend Integer::operator+(int,Integer) */
		friend giv_all_inlined  Integer operator + (const long int l, const Integer& n);
		/** @overload friend Integer::operator+(int,Integer) */
		friend giv_all_inlined  Integer operator + (const long unsigned int l, const Integer& n);
		/** @overload friend Integer::operator+(int,Integer) */
		friend giv_all_inlined  Integer operator + (const Integer& n, const int l);
		/** @overload friend Integer::operator+(int,Integer) */
		friend giv_all_inlined  Integer operator + (const Integer& n, const unsigned int l);


#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
		/** @overload friend Integer::operator+(int,Integer) */
		friend giv_all_inlined 	Integer operator + (const Integer& n, const long long int l);
		/** @overload friend Integer::operator+(int,Integer) */
		friend giv_all_inlined 	Integer operator + (const Integer& n, const long long unsigned int l);
		/** @overload friend Integer::operator+(int,Integer) */
		friend giv_all_inlined 	Integer operator + (const long long int l, const Integer& n);
		/** @overload friend Integer::operator+(int,Integer) */
		friend giv_all_inlined 	Integer operator + (const long long unsigned int l, const Integer& n);
#endif

		//! operator +=.
		//! @param n Integer
		//! @param l to be added up
		friend giv_all_inlined  Integer& operator += (Integer& n, const int l);
		/** @overload friend Integer::operator+=(Integer,int) */
		friend giv_all_inlined  Integer& operator += (Integer& n, const unsigned int l);
		/** @overload friend Integer::operator+=(Integer,int) */
#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
		friend giv_all_inlined 	Integer& operator += (Integer& n, const long long int l);
		/** @overload friend Integer::operator+=(Integer,int) */
		friend giv_all_inlined 	Integer& operator += (Integer& n, const long long unsigned int l);
#endif


		// (FILE gmp++_int_sub.C)

		/*! operator \c -.
		 * @return <code> (*this)-n</code>
		 * @param n as in the formula.
		 */
		giv_all_inlined Integer  operator - (const Integer& n) const;
		/** @overload Integer::operator-(Integer) */
		giv_all_inlined Integer  operator - (const long unsigned int n) const;
		/** @overload Integer::operator-(Integer) */
		giv_all_inlined Integer  operator - (const long int n) const;

		/*! operator \c -= .
		 * @param n as in the formula.
		 * @return <code> (*this) -= n</code>.
		 */
		giv_all_inlined Integer& operator -= (const Integer& n);
		/** @overload Integer::operator-=(Integer) */
		giv_all_inlined Integer& operator -= (const long unsigned int n);
		/** @overload Integer::operator-=(Integer) */
		giv_all_inlined Integer& operator -= (const long int n);
		/** @overload Integer::operator-=(Integer) */
		template<class XXX>
		Integer& operator -=(const XXX& n)
		{
			return this->operator -= ( (Integer)n );
		}

		/*! Opposite.
		 * \return <code>-(*this)</code>.
		 */
		giv_all_inlined Integer  operator -() const;


		//! operator -
		//! @param l,n to be substracted
		friend giv_all_inlined  Integer operator - (const int l, const Integer& n);
		/** @overload friend Integer::operator-(int,Integer) */
		friend giv_all_inlined  Integer operator - (const unsigned int l, const Integer& n);
		/** @overload friend Integer::operator-(int,Integer) */
		friend giv_all_inlined  Integer operator - (const long int l, const Integer& n);
		/** @overload friend Integer::operator-(int,Integer) */
		friend giv_all_inlined  Integer operator - (const long unsigned int l, const Integer& n);
		/** @overload friend Integer::operator-(int,Integer) */
		friend giv_all_inlined  Integer operator - (const Integer& n, const int l);
		/** @overload friend Integer::operator-(int,Integer) */
		friend giv_all_inlined  Integer operator - (const Integer& n, const unsigned int l);

		//! operator -=
		//! @param l,n to be substracted
		friend giv_all_inlined  Integer& operator -= (Integer& n, const int l);
		/** @overload friend Integer::operator-=(Integer,int) */
		friend giv_all_inlined  Integer& operator -= (Integer& n, const unsigned int l);

#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
		/** @overload friend Integer::operator-(int,Integer) */
		friend giv_all_inlined  Integer operator - (const Integer& n, const long long int l);
		/** @overload friend Integer::operator-(int,Integer) */
		friend giv_all_inlined  Integer operator - (const Integer& n, const long long unsigned int l);
		/** @overload friend Integer::operator-(int,Integer) */
		friend giv_all_inlined  Integer operator - (const long long int l, const Integer& n);
		/** @overload friend Integer::operator-(int,Integer) */
		friend giv_all_inlined  Integer operator - (const long long unsigned int l, const Integer& n);

		/** @overload friend Integer::operator-=(Integer,int) */
		friend giv_all_inlined  Integer& operator -= (Integer& n, const long long int l);
		/** @overload friend Integer::operator-=(Integer,int) */
		friend giv_all_inlined  Integer& operator -= (Integer& n, const long long unsigned int l);
#endif

		// (FILE gmp++_int_mul.C)

		/*! operator \c *.
		 * @return <code> (*this)*n</code>
		 * @param n as in the formula.
		 */
		giv_all_inlined Integer  operator * (const Integer& n) const;
		/** @overload Integer::operator*(Integer) */
		giv_all_inlined Integer  operator * (const long unsigned int n) const;
		/** @overload Integer::operator*(Integer) */
		giv_all_inlined Integer  operator * (const long int n) const;

		/*! operator \c *= .
		 * @param n as in the formula.
		 * @return <code> (*this) *= n</code>.
		 */
		giv_all_inlined Integer& operator *= (const Integer& n);
		/** @overload Integer::operator*=(Integer) */
		giv_all_inlined Integer& operator *= (const long unsigned int n);
		/** @overload Integer::operator*=(Integer) */
		giv_all_inlined Integer& operator *= (const long int n);
		/** @overload Integer::operator*=(Integer) */
		template<class XXX>
		Integer& operator *=(const XXX& n) {
			return this->operator *= ( (Integer)n );
		}

		//! operator *
		//! @param l,n to be multpct
		friend giv_all_inlined  Integer operator * (const int l, const Integer& n);
		/** @overload fried Integer::operator*(int,Integer) */
		friend giv_all_inlined  Integer operator * (const unsigned int l, const Integer& n);
		/** @overload fried Integer::operator*(int,Integer) */
		friend giv_all_inlined  Integer operator * (const long int l, const Integer& n);
		/** @overload fried Integer::operator*(int,Integer) */
		friend giv_all_inlined  Integer operator * (const long unsigned int l, const Integer& n);
		/** @overload fried Integer::operator*(int,Integer) */
		friend giv_all_inlined  Integer operator * (const Integer& n, const int l);
		/** @overload fried Integer::operator*(int,Integer) */
		friend giv_all_inlined  Integer operator * (const Integer& n, const unsigned int l);

		//! operator *=
		//! @param l,n to be multpct
		friend giv_all_inlined  Integer& operator *= (Integer& n, const int l);
		/** @overload fried Integer::operator*(Integer,int) */
		friend giv_all_inlined  Integer& operator *= (Integer& n, const unsigned int l);

#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
		/** @overload fried Integer::operator*(int,Integer) */
		friend giv_all_inlined  Integer operator * (const Integer& n, const long long int l);
		/** @overload fried Integer::operator*(int,Integer) */
		friend giv_all_inlined  Integer operator * (const Integer& n, const long long unsigned int l);
		/** @overload fried Integer::operator*(int,Integer) */
		friend giv_all_inlined  Integer operator * (const long long int l, const Integer& n);
		/** @overload fried Integer::operator*(int,Integer) */
		friend giv_all_inlined  Integer operator * (const long long unsigned int l, const Integer& n);

		friend giv_all_inlined  Integer& operator *= (Integer& n, const long long int l);
		/** @overload fried Integer::operator*(Integer,int) */
		friend giv_all_inlined  Integer& operator *= (Integer& n, const long long unsigned int l);
		/** @overload fried Integer::operator*(Integer,int) */
#endif
		///@}

		/*! @name fused add-multiply
		 * @brief Groups a multiplication and an addition/division in a
		 * single function.
		 * This is usually faster than doing the two operations
		 * separately (and preferable to using operators).
		 */
		///@{
		/*! axpy
		 *  <code>res = ax+y</code>.
		 * @param res Integers as in the forumla
		 * @param a Integers as in the forumla
		 * @param x Integers as in the forumla
		 * @param y Integers as in the forumla
		 */
		static giv_all_inlined  Integer& axpy   (Integer& res,
							 const Integer& a,
							 const Integer& x,
							 const Integer& y );
		//! @overload Integer::axpy(Integer,Integer,Integer,Integer)
		static giv_all_inlined  Integer& axpy   (Integer& res,
							 const Integer& a,
							 const long unsigned int x,
							 const Integer& y );

		/*! axpyin (inplace)
		 *  <code>res += ax</code>.
		 * @param res Integers as in the formula.
		 * @param a Integers as in the formula.
		 * @param x Integers as in the formula.
		 */
		static giv_all_inlined  Integer& axpyin   (Integer& res,
							   const Integer& a,
							   const Integer& x);
		//! @overload Integer::axpyin(Integer,Integer,Integer)
		static giv_all_inlined  Integer& axpyin   (Integer& res,
							   const Integer& a,
							   const long unsigned int x);

		/*! maxpy
		 *  <code>res = y - ax</code>.
		 * @param res Integers as in the formula.
		 * @param a Integers as in the formula.
		 * @param x Integers as in the formula.
		 * @param y Integers as in the formula.
		 */
		static giv_all_inlined  Integer& maxpy   (Integer& res,
							  const Integer& a,
							  const Integer& x,
							  const Integer& y );
		//! @overload Integer::maxpy(Integer,Integer,Integer,Integer)
		static giv_all_inlined  Integer& maxpy   (Integer& res,
							  const Integer& a,
							  const long unsigned int x,
							  const Integer& y );

		/*! maxpyin
		 *  <code>res -= ax</code>.
		 * @param res Integers as in the formula.
		 * @param a Integers as in the formula.
		 * @param x Integers as in the formula.
		 */
		static giv_all_inlined  Integer& maxpyin (Integer& res,
							  const Integer& a,
							  const Integer& x );
		//! @overload Integer::maxpyin(Integer,Integer,Integer)
		static giv_all_inlined  Integer& maxpyin (Integer& res,
							  const Integer& a,
							  const long unsigned int x );

		/*! axmy
		 *  <code>res = ax - y</code>.
		 * @param res Integers as in the formula.
		 * @param a Integers as in the formula.
		 * @param x Integers as in the formula.
		 * @param y Integers as in the formula.
		 */
		static giv_all_inlined  Integer& axmy   (Integer& res,
							 const Integer& a,
							 const Integer& x,
							 const Integer& y );
		//! @overload Integer::axmy(Integer,Integer,Integer,Integer)
		static giv_all_inlined  Integer& axmy   (Integer& res,
							 const Integer& a,
							 const long unsigned int x,
							 const Integer & y );

		/*! axmyin (in place)
		 * <code>res = ax - res</code>.
		 * @param res Integers as in the formula.
		 * @param a Integers as in the formula.
		 * @param x Integers as in the formula.
		 */
		static giv_all_inlined  Integer& axmyin   (Integer& res,
							   const Integer& a, const Integer& x);
		//! @overload Integer::axmyin(Integer,Integer,Integer)
		static giv_all_inlined  Integer& axmyin   (Integer& res,
							   const Integer& a, const long unsigned int x);
		///@}

		/*! @name Division/euclidean division/modulo
		 * @brief
		 * The convention for rounding are the following :
		 * -  <code>q = a/b</code>, or equivalent operations with the name \c
		 *   div or \c divin, return \c q  rounded towards \c 0, in the same
		 *   manner as C's '/' (truncated division).
		 * -  <code>r = a % b</code> behaves like C %. The modulo function %
		 *    rounds towards 0 and the sign of the dividend is preserved.  This
		 *    is : \f[ a= b q + r, \text{with } \vert r\vert  < \vert b\vert
		 *    \text{ and } a r \geq 0 \f]
		 * -  <code>r = a mod b</code> or similar functions have the same
		 *    behaviour as GMP \c mpz_mod, that is the remainder is always
		 *    positive (>=0). This is the  division algorithm convention that
		 *    is used (see \c divmod). In a formula : \f[ a= b q + r,
		 *    \text{with } 0 \leq  r  < \vert b\vert \f]
		 *
		 * @warning if <code>q=a/b</code> and <code>r= a % b</code> then <code>
		 * a = b q + r </code> is always true (with in addition <code>0 <= |r|
		 * < |b|</code>).  This is also true for <code>divmod(q,a,b,r)</code>
		 * (and <code>0<=r<|b|</code>).  However, one should not mix the two
		 * conventions and expect equalities <small>(except if a>=0)</small>.
		 */
		// (FILE gmp++_int_div.C)
		///@{
		/*! Division \c q/=d.
		 * @param q quotient
		 * @param d divisor.
		 * @return \c q
		 */
		static giv_all_inlined  Integer& divin (Integer& q, const Integer& d);
		//! @overload Integer::divin(Integer,Integer)
		static giv_all_inlined  Integer& divin (Integer& q, const long int d);
		//! @overload Integer::divin(Integer,Integer)
		static giv_all_inlined  Integer& divin (Integer& q, const long unsigned int d);

		/*! Division \c q=n/d.
		 * @param q quotient
		 * @param n dividand.
		 * @param d divisor
		 * @return \c q
		 */
		static giv_all_inlined  Integer& div   (Integer& q,
							const Integer& n, const Integer& d);
		//! @overload Integer::div(Integer,Integer,Integer)
		static giv_all_inlined  Integer& div   (Integer& q,
							const Integer& n, const long int d);
		//! @overload Integer::div(Integer,Integer,Integer)
		static giv_all_inlined  Integer& div   (Integer& q,
							const Integer& n, const int d);
		//! @overload Integer::div(Integer,Integer,Integer)
		static giv_all_inlined  Integer& div   (Integer& q,
							const Integer& n, const long unsigned int d);

		/*! Division when \c d divides \c n.
		 * @param q exact quotient
		 * @param n dividend
		 * @param d divisor
		 * @warning if quotient is not exact, the result is not predictable.
		 */
		static giv_all_inlined  Integer& divexact (Integer& q,
							   const Integer& n, const Integer& d);
		/*! Division when \c d divides \c n.
		 * @param n dividend
		 * @param d divisor
		 * @return  exact quotient \c n/d
		 * @warning if quotient is not exact, the result is not predictable.
		 */
		static giv_all_inlined  Integer  divexact (const Integer& n, const Integer& d);

		/*! Division operator.
		 * @param d divisor
		 */
		giv_all_inlined Integer  operator /  (const Integer&      d) const;
		//! @overload Integer::operator/(Integer)
		giv_all_inlined Integer  operator /  (const long unsigned int d) const;
		//! @overload Integer::operator/(Integer)
		giv_all_inlined Integer  operator /  (const long int d) const;

		/*! Division operator (inplace).
		 * @param d divisor
		 */
		giv_all_inlined Integer& operator /= (const Integer&      d);
		//! @overload Integer::operator/=(Integer)
		giv_all_inlined Integer& operator /= (const long unsigned int d);
		//! @overload Integer::operator/=(Integer)
		giv_all_inlined Integer& operator /= (const long int d);
		//! @overload Integer::operator/=(Integer)
		template<class XXX>
		Integer& operator /=(const XXX& d) {
			return this->operator /= ( (Integer)d );
		}

		//! operator /
		friend giv_all_inlined  Integer operator / (const int l, const Integer& n);
		//! @overload Integer::operator/(int,Integer)
		friend giv_all_inlined  Integer operator / (const long int l,
							    const Integer& n);
		//! @overload Integer::operator/(int,Integer)
		friend giv_all_inlined  Integer operator / (const Integer& n, const int l);
		//! @overload Integer::operator/(int,Integer)
		friend giv_all_inlined  Integer operator / (const Integer& n, const unsigned int l);

		//! operator /=
		friend giv_all_inlined  Integer& operator /= (Integer& n, const int l);
		//! @overload Integer::operator/=(Integer,int)
		friend giv_all_inlined  Integer& operator /= (Integer& n, const long int l);
		//! @overload Integer::operator/=(Integer,int)
		friend giv_all_inlined  Integer& operator /= (Integer& n, const unsigned int l);

		/*!  Function \c mod (inplace).
		 * \f$ r \gets r \mod n\f$
		 * @param r remainder
		 * @param n modulus
		 */
		static giv_all_inlined  Integer& modin (Integer& r, const Integer& n);
		//! @overload Integer::modin(Integer,Integer)
		static giv_all_inlined  Integer& modin (Integer& r, const long int n);
		//! @overload Integer::modin(Integer,Integer)
		static giv_all_inlined  Integer& modin (Integer& r, const long unsigned int n);
		/*!  Function \c mod.
		 * \f$ r \gets n \mod d\f$
		 * @param r remainder
		 * @param n integer
		 * @param d modulus
		 */
		static giv_all_inlined  Integer& mod   (Integer& r,
							const Integer& n, const Integer& d);
		//! @overload Integer::mod(Integer,Integer,Integer)
		static giv_all_inlined  Integer& mod   (Integer& r,
							const Integer& n, const long int d);
		//! @overload Integer::mod(Integer,Integer,Integer)
		static giv_all_inlined  Integer& mod   (Integer& r,
							const Integer& n, const long unsigned int d);

		/*! Euclidean division.
		 * <code> n = d q + r </code>.
		 * Computes both the quotient and the residue (as in quorem).
		 * @param q as in the formula
		 * @param r as in the formula
		 * @param n as in the formula
		 * @param d as in the formula
		 * @return the quotient \c q
		 */
		static giv_all_inlined  Integer& divmod (Integer& q,
							 Integer& r,
							 const Integer& n,
							 const Integer& d);
		//! @overload Integer::divmod(Integer,Integer,Integer,Integer)
		static giv_all_inlined  Integer& divmod (Integer& q,
							 long int & r,
							 const Integer& n,
							 const long int d);
		//! @overload Integer::divmod(Integer,Integer,Integer,Integer)
		static giv_all_inlined  Integer& divmod (Integer& q,
							 long unsigned int & r,
							 const Integer& n,
							 const long unsigned int d);

		/*! rounding functions.
		 * these are the same as the STL ones, except for the signature.
		 * @param res the result
		 * @param n the numerator
		 * @param d the demominator
		 */
		//! @details same as std::ceil (n/d)
		static giv_all_inlined  Integer&   ceil (Integer & res,
							 const Integer &n,
							 const Integer & d);
		//! @details same as std::floor(n/d)
		static giv_all_inlined  Integer&   floor(Integer & res,
							 const Integer &n,
							 const Integer & d);
		//! @details same as std::trunc(n/d)
		static giv_all_inlined  Integer&   trunc(Integer & res,
							 const Integer &n,
							 const Integer & d);

		/*! rounding functions.
		 * these are the same as the STL ones, except for the signature.
		 * @param n the numerator
		 * @param d the demominator
		 * @return n/d rounded.
		 */
		//! @details same as std::ceil (n/d)
		static giv_all_inlined  Integer    ceil (const Integer &n,
							 const Integer & d);
		//! @details same as std::floor(n/d)
		static giv_all_inlined  Integer    floor(const Integer &n,
							 const Integer & d);
		//! @details same as std::trunc(n/d)
		static giv_all_inlined  Integer    trunc(const Integer &n,
							 const Integer & d);

		// (FILE gmp++_int_mod.C)
		/*! Modulo operator.
		 * @param n modulus
		 * @return remainder <code> (*this) mod n</code>
		 */
		giv_all_inlined Integer  operator % (const Integer& n) const;
		//! @overload Integer::operator%(Integer);
		giv_all_inlined long     operator % (const long unsigned int n) const;
		//! @overload Integer::operator%(Integer);
		giv_all_inlined long     operator % (const long int n) const;
		//! @overload Integer::operator%(Integer);
		giv_all_inlined double   operator % (const double n) const;
		//! @overload Integer::operator%(Integer);
		short    operator % (const unsigned short n) const
		{
			return (short) ( this->operator % ( (long unsigned)n ) );
		}
		//! @overload Integer::operator%(Integer);
		template<class XXX>
		XXX      operator %(const XXX& n) const
		{
			return (XXX)this->operator % ( Integer(n) );
		}

		/** operator %
		 * @param l
		 * @param n
		 * @return n%l
		 */
		friend giv_all_inlined  Integer operator % (const int l, const Integer& n);
		//! @overload Integer::operator%(int,Integer);
		friend giv_all_inlined  Integer operator % (const long int l, const Integer& n);
		//! @overload Integer::operator%(int,Integer);
		friend giv_all_inlined  Integer operator % (const Integer& n, const int l);
		//! @overload Integer::operator%(int,Integer);
		friend giv_all_inlined  Integer operator % (const Integer& n, const unsigned int l);
		//! @overload Integer::operator%(int,Integer);

		friend giv_all_inlined  Integer& operator %= (Integer& n, const int l);
		//! @overload Integer::operator%(int,Integer);
		friend giv_all_inlined  Integer& operator %= (Integer& n, const unsigned int l);


		/*! Modulo operator (inplace).
		 * @param n modulus
		 * @return remainder <code> (*this) <- (*this) mod n</code>
		 */
		giv_all_inlined Integer&  operator %= (const Integer& n);
		//! @overload Integer::operator%=(Integer);
		giv_all_inlined Integer&  operator %= (const long unsigned int n);
		//! @overload Integer::operator%=(Integer);
		giv_all_inlined Integer&  operator %= (const long int n);
#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
		//! @overload Integer::operator%=(Integer);
		Integer&  operator %= (const long long int n)
		{
			return *this %= (Integer)n;
		}
		//! @overload Integer::operator%=(Integer);
		Integer&  operator %= (const long long unsigned int n)
		{
			return *this %= (Integer)n;
		}
		//! @overload Integer::operator%=(Integer);
		giv_all_inlined long long operator % (const long long int n) const;
		//! @overload Integer::operator%=(Integer);
		giv_all_inlined long long unsigned operator % (const long long unsigned int n) const;
#endif
		//! @overload Integer::operator%=(Integer);
		template<class XXX>
		Integer& operator  %=(const XXX& n) {
			return this->operator %= ( (Integer)n );
		}
		///@}






		// (FILE gmp++_int_gcd.C)
		//------------------------------------- Arithmetic functions
		/*! @name Arithmetic functions
		 * @{
		 */
		/** gcd.
		 * @param a,b integers
		 * @return gcd(a,b)
		 */
		friend giv_all_inlined  Integer gcd (const Integer& a, const Integer& b);
		//! \overload Integer::gcd(const Integer&, const Integer&)
		friend giv_all_inlined  Integer gcd ( Integer& u, Integer& v,
						      const Integer& a, const Integer& b);
		//! \overload Integer::gcd(const Integer&, const Integer&)
		friend giv_all_inlined  Integer& gcd (Integer& g, const Integer& a, const Integer& b);
		//! \overload Integer::gcd(const Integer&, const Integer&)
		friend giv_all_inlined  Integer& gcd (Integer& g, Integer& u, Integer& v,
						      const Integer& a, const Integer& b);

		//! Inverse.
		//! Compute the inverse u = a/b.
		/*! @param u
		 * @param a
		 * @param b
		 */
		friend giv_all_inlined  Integer& inv (Integer& u, const Integer& a, const Integer& b);

		//! Compute the inverse inplace  u = u/b.
		/*! @param u
		 * @param b
		 */
		friend giv_all_inlined  Integer& invin (Integer& u, const Integer& b);

		//! pp
		//! @param P,Q params
		friend giv_all_inlined  Integer pp( const Integer& P, const Integer& Q );

		//! lcm
		//! @param g,a,b
		//! @return g=lcm(a,b)
		friend giv_all_inlined  Integer& lcm (Integer& g, const Integer& a, const Integer& b);
		//! lcm
		//! @param a,b
		friend giv_all_inlined  Integer lcm (const Integer& a, const Integer& b);

		// (FILE gmp++_int_pow.C)
		//! pow. return \f$n^l\f$
		//!@param Res,n,l
		friend giv_all_inlined  Integer& pow(Integer& Res, const Integer& n, const long int l);
		//! @overload Integer::pow(Integer,Integer,long int)
		friend giv_all_inlined  Integer& pow(Integer& Res, const long unsigned int n, const long unsigned int l);
		//! @overload Integer::pow(Integer,Integer,long int)
		friend giv_all_inlined  Integer& pow(Integer& Res, const Integer& n, const long unsigned int l);
		//! @overload Integer::pow(Integer,Integer,long int)
		friend giv_all_inlined  Integer& pow(Integer& Res, const Integer& n, const int l)
		{
			return pow(Res, n, (long)l );
		}
		//! @overload Integer::pow(Integer,Integer,long int)
		friend giv_all_inlined  Integer& pow(Integer& Res, const Integer& n, const unsigned int l)
		{
			return pow(Res, n, (long unsigned)l );
		}

		//! pow. return \f$n^l\f$
		//!@param n,l
		friend giv_all_inlined  Integer pow(const Integer& n, const long int l);
		//! @overload Integer::pow(Integer,long int)
		friend giv_all_inlined  Integer pow(const Integer& n, const long unsigned int l);
		//! @overload Integer::pow(Integer,long int)
		friend giv_all_inlined  Integer pow(const Integer& n, const int l)
		{
			return pow(n, (long)l );
		}
		//! @overload Integer::pow(Integer,long int)
		friend giv_all_inlined  Integer pow(const Integer& n, const unsigned int l)
		{
			return pow(n, (long unsigned)l );
		}
		///@}

		//! modular pow. return \f$n^e \mod  m\f$.
		friend giv_all_inlined  Integer& powmod(Integer& Res, const Integer& n, const long unsigned int e, const Integer& m);
		//! @overload Integer::powmod(Integer,Integer,long unsigned int,Integer)
		friend giv_all_inlined  Integer& powmod(Integer& Res, const Integer& n, const long int e, const Integer& m);
		//! @overload Integer::powmod(Integer,Integer,long unsigned int,Integer)
		friend giv_all_inlined  Integer& powmod(Integer& Res, const Integer& n, const unsigned int e, const Integer& m)
		{
			return powmod(Res, n, (long unsigned)e, m);
		}
		//! @overload Integer::powmod(Integer,Integer,long unsigned int,Integer)
		friend giv_all_inlined  Integer& powmod(Integer& Res, const Integer& n, const int e, const Integer& m)
		{
			return powmod(Res, n, (long)e, m);
		}
		//! @overload Integer::powmod(Integer,Integer,long unsigned int,Integer)
		friend giv_all_inlined  Integer& powmod(Integer& Res, const Integer& n, const Integer& e, const Integer& m);


		//! modular pow. return \f$n^e \mod  m\f$.
		//! @param n,e,m
		friend giv_all_inlined  Integer powmod(const Integer& n, const long unsigned int e, const Integer& m);
		//! @overload Integer::powmod(Integer,long unsigned int,Integer)
		friend giv_all_inlined  Integer powmod(const Integer& n, const long int e, const Integer& m);
		//! @overload Integer::powmod(Integer,long unsigned int,Integer)
		friend giv_all_inlined  Integer powmod(const Integer& n, const unsigned int e, const Integer& m)
		{
			return powmod(n, (long unsigned)e, m);
		}
		//! @overload Integer::powmod(Integer,long unsigned int,Integer)
		friend giv_all_inlined  Integer powmod(const Integer& n, const int e, const Integer& m)
		{
			return powmod(n, (long)e, m);
		}
		//! @overload Integer::powmod(Integer,long unsigned int,Integer)
		friend giv_all_inlined  Integer powmod(const Integer& n, const Integer& e, const Integer& m);

		// (FILE gmp++_int_misc.C)
		//! fact
		//! @param l
		friend giv_all_inlined  Integer fact ( long unsigned int l);

		//! (square) roots
		//! @param p
		friend giv_all_inlined  Integer sqrt(const Integer& p);
		//! (square) roots
		//! @param r,p
		friend giv_all_inlined  Integer& sqrt(Integer& r, const Integer& p);

		//! (square) roots
		//! @param p,rem
		friend giv_all_inlined  Integer sqrtrem(const Integer& p, Integer& rem);
		//! (square) roots
		//! @param r,p,rem
		friend giv_all_inlined  Integer& sqrtrem(Integer& r, const Integer& p, Integer& rem);

		//! (square) roots
		//! @param q,a,n
		friend giv_all_inlined  bool root(Integer& q,
						  const Integer&a, unsigned int n);

		//! logs
		//! @param a,p
		friend giv_all_inlined  long logp(const Integer& a, const Integer& p) ;
		//! logs
		//! @param a
		friend giv_all_inlined  double logtwo(const Integer& a) ;
		//! logs
		//! @param a
		friend giv_all_inlined  double naturallog(const Integer& a) ;
		///@}

		//-----------------------------------------Miscellaneous
		/*! @name Miscellaneous.
		*/
		///@{
		//! swap
		//! @param a,b
		friend giv_all_inlined  void swap(Integer& a, Integer&b);

		//! sign
		//! @param a
 		friend inline int sign   (const Integer& a)
		{
			return a.priv_sign();
		}

		//! sign
		int sign() const // BB is this usefull ?
		{
			return priv_sign();
		} // but figure out the friend sign()

		//! private sign
		int priv_sign() const // BB no longer protected.
		{
			return mpz_sgn( (mpz_srcptr)&gmp_rep );
		}
		///@}


		//! @name representation
		///@{
		//! get representation
		mpz_ptr get_mpz()
		{
			return (mpz_ptr)&gmp_rep;
		}

		//! get representation (constant)
		mpz_srcptr get_mpz_const() const
		{
			return (mpz_srcptr)&gmp_rep;
		}

		/*! returns the number of bytes used to store *this          */
		//! @param a
		friend giv_all_inlined  long unsigned int length (const Integer& a);

		/*! returns the number of machine words used to store *this          */
		giv_all_inlined size_t size() const;

		/*! returns <code>ceil(log_BASE(*this))</code>.        */
		giv_all_inlined size_t size_in_base(int B) const;

		/*! returns <code>ceil(log_2(*this)) </code>.        */
		giv_all_inlined size_t bitsize() const;

		/*! return the i-th word of the integer. Word 0 is lowest word.*/
		giv_all_inlined long unsigned operator[](size_t i) const;



		// (FILE gmp++_int_misc.C)

		//! perfect power
		friend giv_all_inlined  int isperfectpower  (const Integer& n );

		//! absolute value
		friend giv_all_inlined  Integer abs(const Integer& n);

		//! @name primes
		///@{
		friend giv_all_inlined  Integer& prevprime(Integer&, const Integer& p);
		friend giv_all_inlined  Integer& nextprime(Integer&, const Integer& p);
		friend giv_all_inlined  int probab_prime(const Integer& p);
		friend giv_all_inlined  int probab_prime(const Integer& p, int r);
		friend giv_all_inlined  int jacobi(const Integer& u, const Integer& v) ;
		friend giv_all_inlined  int legendre(const Integer& u, const Integer& v) ;
		///@}

		//! @name Increment/Decrement operators
		///@{
		Integer& operator++()
		{ // prefix
			return *this+=1UL;
		}
		Integer operator++(int)
		{ // postfix
			Integer tmp = *this ;
			++*this;
			return tmp;
		}
		Integer& operator--()
		{// prefix
			return *this-=1UL;
		}
		Integer operator--(int)
		{// postfix
			Integer tmp = *this ;
			--*this;
			return tmp;
		}
		///@}


		/** @name Cast operators.
		 * Convert an Integer to a basic C++ type.
		 * @warning Cast towards \b unsigned consider only the absolute
		 * value
		 */
		///@{
		operator bool() const
		{
			return *this!=0UL;
		}
		operator short() const
		{
			return (short)(int) *this;
		}
		operator unsigned short() const
		{
			return (unsigned short) (unsigned int) *this;
		}
		operator unsigned char() const
		{
			return (unsigned char) (unsigned int) *this;
		}
		giv_all_inlined operator unsigned int() const ;
		giv_all_inlined operator int() const ;
		operator signed char() const
		{
			return (signed char) (int) *this;
		}
		giv_all_inlined operator long unsigned() const ;
		giv_all_inlined operator long() const ;
#ifndef __GIVARO__DONOTUSE_longlong__
		giv_all_inlined operator long long unsigned() const ;
		giv_all_inlined operator long long() const ;
#endif
		giv_all_inlined operator std::string() const ;
		giv_all_inlined operator float() const ;
		giv_all_inlined operator double() const ;
		giv_all_inlined operator vect_t() const ;
		///@}

		// (FILE gmp++_int_rand.inl)
		/* BB:
		 * if the following functions are NOT static inline, one
		 * can use -Wl,-zmuldefs....
		 */
		//--------------------Random Iterators
		/*! @name Random numbers functions
		*/
		//! Random numbers (no doc)
		///@{
		static inline void seeding(long unsigned int  s);
		static inline void seeding(Integer s);
		static inline void seeding();

#ifdef __GMP_PLUSPLUS__
		static inline gmp_randclass& randstate();
#else
		// static __gmp_randstate_struct initializerandstate();
		// static __gmp_randstate_struct* randstate();
#endif
		static inline bool RandBool()  ;
		/*  random <= */
		template<bool ALWAYSPOSITIVE>
		static inline Integer& random_lessthan (Integer& r, const Integer & m);
		static inline Integer& random_lessthan (Integer& r, const Integer & m) ;
		template<bool ALWAYSPOSITIVE>
		static inline Integer& random_lessthan_2exp (Integer& r, const long unsigned int & m);

		static inline Integer& random_lessthan_2exp (Integer& r, const long unsigned int & m) ;
		template<bool ALWAYSPOSITIVE>
		static inline Integer random_lessthan_2exp (const long unsigned int & m);
		static inline Integer random_lessthan_2exp (const long unsigned int & m) ;
		template<bool ALWAYSPOSITIVE>
		static inline Integer& random_lessthan (Integer& r, const long unsigned int & m) ;

		static inline Integer& random_lessthan (Integer& r, const long unsigned int & m) ;
		template<bool ALWAYSPOSITIVE,class T>
		static inline Integer random_lessthan (const T & m);
		template<class T>
		static inline Integer random_lessthan (const T & m) ;


		/*  random = */
		template<bool ALWAYSPOSITIVE>
		static inline Integer& random_exact_2exp (Integer& r, const long unsigned int & m) ;
		static inline Integer& random_exact_2exp (Integer& r, const long unsigned int & m);


		template<bool ALWAYSPOSITIVE>
		static inline Integer& random_exact (Integer& r, const Integer & s) ;
		static inline Integer& random_exact (Integer& r, const Integer & s) ;

		template<bool ALWAYSPOSITIVE>
		static inline Integer& random_exact (Integer& r, const long unsigned int & m)  ;
		static inline Integer& random_exact (Integer& r, const long unsigned int & m) ;
		template<bool ALWAYSPOSITIVE,class T>
		static inline Integer& random_exact (Integer& r, const T & m) ;
		template<class T>
		static inline Integer& random_exact (Integer& r, const T & m) ;
		template<bool ALWAYSPOSITIVE,class T>
		static inline Integer random_exact (const T & s) ;
		template<class T>
		static inline Integer random_exact (const T & s) ;

		/*  random <.< */
		static inline Integer& random_between (Integer& r, const Integer& m, const Integer&M) ;
		static inline Integer random_between (const Integer& m, const Integer &M) ;
		static inline Integer& random_between_2exp (Integer& r, const long unsigned int& m,
							    const long unsigned int &M) ;
		static inline Integer& random_between (Integer& r, const long unsigned int& m,
						       const long unsigned int &M) ;
		static inline Integer random_between_2exp (const long unsigned int & m,
							   const long unsigned int &M) ;
		static inline Integer random_between (const long unsigned int & m,
						      const long unsigned int &M) ;
		template<class R>
		static inline Integer random_between (const R & m, const R & M) ;
		template<class R>
		static inline Integer & random_between (Integer &r, const R & m, const R & M);


		// useful functions :
		template<bool ALWAYSPOSITIVE,class T>
		static inline Integer& random (Integer& r, const T & m) ;
		template<class T>
		static inline Integer& random (Integer& r, const T & m) ;
		template<bool ALWAYSPOSITIVE,class T>
		static inline Integer random(const T & sz) ;
		template<class T>
		static inline Integer random(const T & sz) ;
		template<bool ALWAYSPOSITIVE>
		static inline Integer random();
		static inline Integer random();
		template<bool ALWAYSPOSITIVE,class T>
		static inline Integer nonzerorandom(const T & sz) ;
		template<bool ALWAYSPOSITIVE,class T>
		static inline Integer& nonzerorandom (Integer& r, const T& size) ;
		template<class T>
		static inline Integer nonzerorandom(const T & sz) ;
		template<class T>
		static inline Integer& nonzerorandom (Integer& r, const T& size) ;
		static inline Integer nonzerorandom() ;
		///@}

		//----------------------------------------------I/O
		/*! @name I/O */
		//! Input/Output of Integers
		///@{

		/** in operator.
		 * @param i input stream
		 * @param n integer to be built
		 */
		friend giv_all_inlined  std::istream& operator >> (std::istream &i, Integer& n);

		/** @brief out operator.
		 * @param o output stream
		 * @param n integer to be printed
		 */
		friend giv_all_inlined  std::ostream& operator << (std::ostream &o, const Integer& n);

		//! nodoc
		//! @param o output
		//! @param n integer
		friend giv_all_inlined  std::ostream& absOutput (std::ostream &o, const Integer& n);

		/// nodoc
		/*! @param  x x
		 * @param count x
		 * @param order x
		 * @param size x
		 * @param endian x
		 * @param nails x
		 * @param op x
		 */
		friend giv_all_inlined  void importWords(Integer& x , size_t count, int order,
							 int size, int endian, size_t nails,
							 const void* op);
		/** print integer.
		 * @param o output stream.
		 */
		giv_all_inlined std::ostream& print( std::ostream& o ) const;
		///@}


	protected:

		typedef __mpz_struct Rep; //!< @internal rep type

		Rep gmp_rep;//!< @internal rep

		//mpz_ptr get_mpz()
		//{ return (mpz_ptr)&gmp_rep; }

		//!@internal get representation.
		const Rep* get_rep() const
		{
			return &gmp_rep;
		}


	}; //----------------------------------------------- End of Class Integer

	//! generic sign
	template<class T>
	static inline int sign (const T a)
	{
		// std::cout << ((a>0)-(a<0)) << std::endl;
		return (a>0)-(a<0);
		// if (a == 0) return 0 ;
		// return (a<0)?-1:1 ;
	}


} // Givaro

// only template code is inlined
#ifdef __GIVARO_INLINE_ALL
#include "gmp++/gmp++_int.C"
#endif
#include "gmp++/gmp++_int_rand.inl"

#endif // __GIVARO_GMPplusplus_integer_H

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen
