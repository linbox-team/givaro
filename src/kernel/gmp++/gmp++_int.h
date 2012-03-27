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
	Integer& 	pow(Integer& Res, const Integer& n, const long l);
	Integer& 	pow(Integer& Res, const unsigned long n, const unsigned long l);
	Integer& 	pow(Integer& Res, const Integer& n, const unsigned long l);
	Integer& 	pow(Integer& Res, const Integer& n, const int l) ;
	Integer& 	pow(Integer& Res, const Integer& n, const unsigned int l) ;
	Integer 	pow(const Integer& n, const long l);
	Integer 	pow(const Integer& n, const unsigned long l);
	Integer 	pow(const Integer& n, const int l) ;
	Integer 	pow(const Integer& n, const unsigned int l);
	Integer& 	powmod(Integer& Res, const Integer& n, const unsigned long e, const Integer& m);
	Integer& 	powmod(Integer& Res, const Integer& n, const long e, const Integer& m);
	Integer& 	powmod(Integer& Res, const Integer& n, const unsigned int e, const Integer& m) ;
	Integer& 	powmod(Integer& Res, const Integer& n, const int e, const Integer& m)  ;
	Integer& 	powmod(Integer& Res, const Integer& n, const Integer& e, const Integer& m);
	Integer 	powmod(const Integer& n, const unsigned long e, const Integer& m);
	Integer 	powmod(const Integer& n, const long e, const Integer& m);
	Integer 	powmod(const Integer& n, const unsigned int e, const Integer& m) ;
	Integer 	powmod(const Integer& n, const int e, const Integer& m) ;
	Integer 	powmod(const Integer& n, const Integer& e, const Integer& m);
	// (FILE inline)
	int 		sign   (const Integer& a);
	// (FILE gmp++_int_compare.C)
	int 		isZero (const Integer& a);
	int 		compare(const Integer& a, const Integer& b);
	int 		absCompare(const Integer& a, const Integer& b);
	int 		absCompare(const Integer& a, const double d);
	int 		absCompare(const Integer& a, const float d);
	int 		absCompare(const Integer& a, const unsigned long u);
	int 		absCompare(const Integer& a, const unsigned u);
	int 		absCompare(const Integer& a, const long int u);
	int 		absCompare(const Integer& a, const int u);
	template<class T>
	int             absCompare( const T a,  const Integer & b)
	{
		return absCompare(b,a);
	}
	// template<class A, class B>
	// bool isleq(const A&a, const B&b);
	int 		nonZero (const Integer& a);
	int 		isOne  (const Integer& a);
	// (FILE gmp++_int_misc.C)
	Integer 	fact ( unsigned long l);
	Integer 	sqrt(const Integer& p);
	Integer 	sqrtrem(const Integer& p, Integer& rem);
	Integer& 	sqrt(Integer& r, const Integer& p);
	Integer& 	sqrtrem(Integer& r, const Integer& p, Integer& rem);
	bool 		root(Integer& q, const Integer&, unsigned int n);
	long 		logp(const Integer& a, const Integer& p) ;
	double 		logtwo(const Integer& a) ;
	void 		swap(Integer& , Integer&);
	double 		naturallog(const Integer& a) ;
	int 		isperfectpower  (const Integer& );
	Integer 	abs(const Integer& n);
	Integer& 	prevprime(Integer&, const Integer& p);
	Integer& 	nextprime(Integer&, const Integer& p);
	int 		probab_prime(const Integer& p);
	int 		probab_prime(const Integer& p, int r);
	int 		jacobi(const Integer& u, const Integer& v) ;
	int 		legendre(const Integer& u, const Integer& v) ;
	unsigned long 	length (const Integer& a);
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
		/*! @name Constructor/Destructors
		*/
		//@{
		giv_all_inlined Integer( const std::vector<mp_limb_t>& vect_t );
		giv_all_inlined Integer(int n = 0);
		giv_all_inlined Integer(long n);
		giv_all_inlined Integer(unsigned char n);
		giv_all_inlined Integer(unsigned int n);
		giv_all_inlined Integer(unsigned long n);
#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
		giv_all_inlined Integer(long long n);
		giv_all_inlined Integer(unsigned long long n);
#endif
		//! Creates a new Integer.
		/*! @param sz size
		 * @param d array
		 */
		giv_all_inlined Integer(unsigned long* d, long sz);

		giv_all_inlined Integer(double d);
		giv_all_inlined Integer(const char *s);
		giv_all_inlined Integer(const Integer& n);
		giv_all_inlined ~Integer();
		//@}

		//------------------------------------- predefined null and one
		//! zero
		static const Integer zero ;
		//! one
		static const Integer one ;
		static const Integer mOne ;


		// -- Assignment and copy operators
		/*! @name Assignment and copy operators
		*/
		//@{
		giv_all_inlined Integer& operator = (const Integer& n);
		giv_all_inlined Integer& logcpy(const Integer& n);
		giv_all_inlined Integer& copy(const Integer& n);
		//@}

		//------------------Equalities and inequalities between integers and longs
		// (FILE gmp++_int_compare.C)
		/*! @name (in)equality
		*/
		//@{

		giv_all_inlined friend int compare(const Integer& a, const Integer& b);
		giv_all_inlined friend int absCompare(const Integer& a, const Integer& b);
		giv_all_inlined friend int absCompare(const Integer& a, const double d);
		giv_all_inlined friend int absCompare(const Integer& a, const float d);
		giv_all_inlined friend int absCompare(const Integer& a, const unsigned long u);
		giv_all_inlined friend int absCompare(const Integer& a, const unsigned u);
		giv_all_inlined friend int absCompare(const Integer& a, const long int u);
		giv_all_inlined friend int absCompare(const Integer& a, const int u);
		template<class T>
		giv_all_inlined friend int absCompare( const T a,  const Integer & b) ;


		giv_all_inlined int operator >= (const Integer & l) const;
		giv_all_inlined int operator >= (const int l) const;
		giv_all_inlined int operator >= (const long l) const;
		giv_all_inlined int operator >= (const unsigned long l) const;
		giv_all_inlined int operator >= (const unsigned l) const;
		giv_all_inlined int operator >= (const double l) const;
		giv_all_inlined int operator >= (const float l) const;
		giv_all_inlined friend int operator >= (float l, const Integer& n);
		giv_all_inlined friend int operator >= (double l, const Integer& n);
		giv_all_inlined friend int operator >= (int l, const Integer& n);
		giv_all_inlined friend int operator >= (long l, const Integer& n);
		giv_all_inlined friend int operator >= (unsigned l, const Integer& n);
		giv_all_inlined friend int operator >= (unsigned long l, const Integer& n);


		giv_all_inlined int operator <= (const Integer & l) const;
		giv_all_inlined int operator <= (const int l) const;
		giv_all_inlined int operator <= (const long l) const;
		giv_all_inlined int operator <= (const unsigned long l) const;
		giv_all_inlined int operator <= (const unsigned l) const;
		giv_all_inlined int operator <= (const double l) const;
		giv_all_inlined int operator <= (const float l) const;
		giv_all_inlined friend int operator <= (float l, const Integer& n);
		giv_all_inlined friend int operator <= (double l, const Integer& n);
		giv_all_inlined friend int operator <= (int l, const Integer& n);
		giv_all_inlined friend int operator <= (long l, const Integer& n);
		giv_all_inlined friend int operator <= (unsigned l, const Integer& n);
		giv_all_inlined friend int operator <= (unsigned long l, const Integer& n);


		giv_all_inlined int operator != (const Integer & l) const;
		giv_all_inlined int operator != (const int l) const;
		giv_all_inlined int operator != (const long l) const;
		giv_all_inlined int operator != (const unsigned long l) const;
		giv_all_inlined int operator != (const unsigned l) const;
		giv_all_inlined int operator != (const double l) const;
		giv_all_inlined int operator != (const float l) const;
		giv_all_inlined friend int operator != (float l, const Integer& n);
		giv_all_inlined friend int operator != (double l, const Integer& n);
		giv_all_inlined friend int operator != (int l, const Integer& n);
		giv_all_inlined friend int operator != (long l, const Integer& n);
		giv_all_inlined friend int operator != (unsigned l, const Integer& n);
		giv_all_inlined friend int operator != (unsigned long l, const Integer& n);


		giv_all_inlined int operator == (const Integer & l) const;
		giv_all_inlined int operator == (const int l) const;
		giv_all_inlined int operator == (const long l) const;
		giv_all_inlined int operator == (const unsigned long l) const;
		giv_all_inlined int operator == (const unsigned l) const;
		giv_all_inlined int operator == (const double l) const;
		giv_all_inlined int operator == (const float l) const;
		giv_all_inlined friend int operator == (float l, const Integer& n);
		giv_all_inlined friend int operator == (double l, const Integer& n);
		giv_all_inlined friend int operator == (int l, const Integer& n);
		giv_all_inlined friend int operator == (long l, const Integer& n);
		giv_all_inlined friend int operator == (unsigned l, const Integer& n);
		giv_all_inlined friend int operator == (unsigned long l, const Integer& n);


		giv_all_inlined int operator > (const Integer & l) const;
		giv_all_inlined int operator > (const int l) const;
		giv_all_inlined int operator > (const long l) const;
		giv_all_inlined int operator > (const unsigned long l) const;
		giv_all_inlined int operator > (const unsigned l) const;
		giv_all_inlined int operator > (const double l) const;
		giv_all_inlined int operator > (const float l) const;
		giv_all_inlined friend int operator > (float l, const Integer& n);
		giv_all_inlined friend int operator > (double l, const Integer& n);
		giv_all_inlined friend int operator > (int l, const Integer& n);
		giv_all_inlined friend int operator > (long l, const Integer& n);
		giv_all_inlined friend int operator > (unsigned l, const Integer& n);
		giv_all_inlined friend int operator > (unsigned long l, const Integer& n);

		giv_all_inlined int operator < (const Integer & l) const;
		giv_all_inlined int operator < (const int l) const;
		giv_all_inlined int operator < (const long l) const;
		giv_all_inlined int operator < (const unsigned long l) const;
		giv_all_inlined int operator < (const unsigned l) const;
		giv_all_inlined int operator < (const double l) const;
		giv_all_inlined int operator < (const float l) const;
		giv_all_inlined friend int operator < (float l, const Integer& n);
		giv_all_inlined friend int operator < (double l, const Integer& n);
		giv_all_inlined friend int operator < (int l, const Integer& n);
		giv_all_inlined friend int operator < (long l, const Integer& n);
		giv_all_inlined friend int operator < (unsigned l, const Integer& n);
		giv_all_inlined friend int operator < (unsigned long l, const Integer& n);





		//@}

		//------------------ Bit logic
		// (FILE gmp++_int_misc.C)
		/*!@name Bit logic
		*/
		//@{
		giv_all_inlined Integer operator^ (const Integer&) const;   // XOR
		giv_all_inlined Integer operator| (const Integer&) const;   // OR
		giv_all_inlined Integer operator& (const Integer&) const;   // AND
		giv_all_inlined unsigned long operator^ (const unsigned long& a) const;
		giv_all_inlined unsigned long operator| (const unsigned long& a) const;
		giv_all_inlined unsigned long operator& (const unsigned long& a) const;
		giv_all_inlined Integer operator ~ () const;   // 1 complement
		giv_all_inlined Integer& operator^= (const Integer&);   // XOR
		giv_all_inlined Integer& operator|= (const Integer&);   // OR
		giv_all_inlined Integer& operator&= (const Integer&);   // AND
		giv_all_inlined Integer operator<< (int l) const; // lshift
		giv_all_inlined Integer operator>> (int l) const; // rshift
		giv_all_inlined Integer operator<< (long l) const; // lshift
		giv_all_inlined Integer operator>> (long l) const; // rshift
		giv_all_inlined Integer operator<< (unsigned int l) const; // lshift
		giv_all_inlined Integer operator>> (unsigned int l) const; // rshift
		giv_all_inlined Integer operator<< (unsigned long l) const; // lshift
		giv_all_inlined Integer operator>> (unsigned long l) const; // rshift
		giv_all_inlined Integer& operator<<= (int l) ; // lshift
		giv_all_inlined Integer& operator>>= (int l) ; // rshift
		giv_all_inlined Integer& operator<<= (long l) ; // lshift
		giv_all_inlined Integer& operator>>= (long l) ; // rshift
		giv_all_inlined Integer& operator<<= (unsigned int l) ; // lshift
		giv_all_inlined Integer& operator>>= (unsigned int l) ; // rshift
		giv_all_inlined Integer& operator<<= (unsigned long l) ; // lshift
		giv_all_inlined Integer& operator>>= (unsigned long l) ; // rshift
		//@}

		// - Methods
		/*! @name Addition, substraction, multiplication
		*/
		//@{
		// (FILE gmp++_int_add.C)
		/*!  Addition (inplace)
		 * <code>res+=n</code>.
		 * @param res as in the formula
		 * @param n as in the formula
		 */
		static giv_all_inlined  Integer& addin (Integer& res, const Integer& n);
		static giv_all_inlined  Integer& addin (Integer& res, const long n);
		static giv_all_inlined  Integer& addin (Integer& res, const unsigned long n);
		/*!  Addition
		 * <code>res=n1+n2</code>.
		 * @param res as in the formula
		 * @param n1 as in the formula
		 * @param n2 as in the formula
		 */
		static giv_all_inlined  Integer& add   (Integer& res, const Integer& n1, const Integer& n2);
		static giv_all_inlined  Integer& add   (Integer& res, const Integer& n1, const long n2);
		static giv_all_inlined  Integer& add   (Integer& res, const Integer& n1, const unsigned long n2);

		// (FILE gmp++_int_sub.C)
		/*!  Substraction (inplace)
		 * <code>res-=n</code>.
		 * @param res as in the formula
		 * @param n as in the formula
		 */
		static giv_all_inlined  Integer& subin (Integer& res, const Integer& n);
		static giv_all_inlined  Integer& subin (Integer& res, const long n);
		static giv_all_inlined  Integer& subin (Integer& res, const unsigned long n);
		/*!  Substraction
		 * <code>res=n1-n2</code>.
		 * @param res as in the formula
		 * @param n1 as in the formula
		 * @param n2 as in the formula
		 */
		static giv_all_inlined  Integer& sub   (Integer& res, const Integer& n1, const Integer& n2);
		static giv_all_inlined  Integer& sub   (Integer& res, const Integer& n1, const long n2);
		static giv_all_inlined  Integer& sub   (Integer& res, const Integer& n1, const unsigned long n2);
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
		static giv_all_inlined  Integer& mulin (Integer& res, const long n);
		static giv_all_inlined  Integer& mulin (Integer& res, const unsigned long n);
		/*! Multiplication
		 * <code>res=n1*n2</code>.
		 * @param res as in the formula
		 * @param n1 as in the formula
		 * @param n2 as in the formula
		 */
		static giv_all_inlined  Integer& mul   (Integer& res, const Integer& n1, const Integer& n2);
		static giv_all_inlined  Integer& mul   (Integer& res, const Integer& n1, const long n2);
		static giv_all_inlined  Integer& mul   (Integer& res, const Integer& n1, const unsigned long n2);

		//----------------Elementary arithmetic between Integers & longs
		// (FILE gmp++_int_add.C)
		/*! operator \c +.
		 * @return <code> (*this)+n</code>
		 * @param n as in the formula.
		 */
		giv_all_inlined Integer  operator + (const Integer& n) const;
		giv_all_inlined Integer  operator + (const unsigned long n) const;
		giv_all_inlined Integer  operator + (const long n) const;
		/*! operator \c += .
		 * @param n asfriend In the formula.
		 * @return <code> (*this) += n</code>.
		 */
		giv_all_inlined Integer& operator += (const Integer& n);
		giv_all_inlined Integer& operator += (const unsigned long n);
		giv_all_inlined Integer& operator += (const long n);
		template<class XXX>
		Integer& operator +=(const XXX& n) {
			return this->operator += ( (Integer)n );
		}

		friend giv_all_inlined  Integer operator + (const int l, const Integer& n);
		friend giv_all_inlined  Integer operator + (const unsigned int l, const Integer& n);
		friend giv_all_inlined  Integer operator + (const long l, const Integer& n);
		friend giv_all_inlined  Integer operator + (const unsigned long l, const Integer& n);
		friend giv_all_inlined  Integer operator + (const Integer& n, const int l);
		friend giv_all_inlined  Integer operator + (const Integer& n, const unsigned int l);

		friend giv_all_inlined  Integer& operator += (Integer& n, const int l);
		friend giv_all_inlined  Integer& operator += (Integer& n, const unsigned int l);

#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
		friend giv_all_inlined 	Integer operator + (const Integer& n, const long long l);
		friend giv_all_inlined 	Integer operator + (const Integer& n, const unsigned long long l);
		friend giv_all_inlined 	Integer operator + (const long long l, const Integer& n);
		friend giv_all_inlined 	Integer operator + (const unsigned long long l, const Integer& n);
		friend giv_all_inlined 	Integer& operator += (Integer& n, const long long l);
		friend giv_all_inlined 	Integer& operator += (Integer& n, const unsigned long long l);
#endif


		// (FILE gmp++_int_sub.C)

		/*! operator \c -.
		 * @return <code> (*this)-n</code>
		 * @param n as in the formula.
		 */
		giv_all_inlined Integer  operator - (const Integer& n) const;
		giv_all_inlined Integer  operator - (const unsigned long n) const;
		giv_all_inlined Integer  operator - (const long n) const;
		/*! operator \c -= .
		 * @param n as in the formula.
		 * @return <code> (*this) -= n</code>.
		 */
		giv_all_inlined Integer& operator -= (const Integer& n);
		giv_all_inlined Integer& operator -= (const unsigned long n);
		giv_all_inlined Integer& operator -= (const long n);
		template<class XXX>
		Integer& operator -=(const XXX& n)
		{
			return this->operator -= ( (Integer)n );
		}

		/*! Opposite.
		 * \return <code>-(*this)</code>.
		 */
		giv_all_inlined Integer  operator -() const;

		friend giv_all_inlined  Integer operator - (const int l, const Integer& n);
		friend giv_all_inlined  Integer operator - (const unsigned int l, const Integer& n);
		friend giv_all_inlined  Integer operator - (const long l, const Integer& n);
		friend giv_all_inlined  Integer operator - (const unsigned long l, const Integer& n);
		friend giv_all_inlined  Integer operator - (const Integer& n, const int l);
		friend giv_all_inlined  Integer operator - (const Integer& n, const unsigned int l);

		friend giv_all_inlined  Integer& operator -= (Integer& n, const int l);
		friend giv_all_inlined  Integer& operator -= (Integer& n, const unsigned int l);

#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
		friend giv_all_inlined  Integer operator - (const Integer& n, const long long l);
		friend giv_all_inlined  Integer operator - (const Integer& n, const unsigned long long l);
		friend giv_all_inlined  Integer operator - (const long long l, const Integer& n);
		friend giv_all_inlined  Integer operator - (const unsigned long long l, const Integer& n);

		friend giv_all_inlined  Integer& operator -= (Integer& n, const long long l);
		friend giv_all_inlined  Integer& operator -= (Integer& n, const unsigned long long l);
#endif

		// (FILE gmp++_int_mul.C)

		/*! operator \c *.
		 * @return <code> (*this)*n</code>
		 * @param n as in the formula.
		 */
		giv_all_inlined Integer  operator * (const Integer& n) const;
		giv_all_inlined Integer  operator * (const unsigned long n) const;
		giv_all_inlined Integer  operator * (const long n) const;
		/*! operator \c *= .
		 * @param n as in the formula.
		 * @return <code> (*this) *= n</code>.
		 */
		giv_all_inlined Integer& operator *= (const Integer& n);
		giv_all_inlined Integer& operator *= (const unsigned long n);
		giv_all_inlined Integer& operator *= (const long n);
		template<class XXX>
		Integer& operator *=(const XXX& n) {
			return this->operator *= ( (Integer)n );
		}

		// -- operator *
		friend giv_all_inlined  Integer operator * (const int l, const Integer& n);
		friend giv_all_inlined  Integer operator * (const unsigned int l, const Integer& n);
		friend giv_all_inlined  Integer operator * (const long l, const Integer& n);
		friend giv_all_inlined  Integer operator * (const unsigned long l, const Integer& n);
		friend giv_all_inlined  Integer operator * (const Integer& n, const int l);
		friend giv_all_inlined  Integer operator * (const Integer& n, const unsigned int l);

		friend giv_all_inlined  Integer& operator *= (Integer& n, const int l);
		friend giv_all_inlined  Integer& operator *= (Integer& n, const unsigned int l);

#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
		friend giv_all_inlined  Integer operator * (const Integer& n, const long long l);
		friend giv_all_inlined  Integer operator * (const Integer& n, const unsigned long long l);
		friend giv_all_inlined  Integer operator * (const long long l, const Integer& n);
		friend giv_all_inlined  Integer operator * (const unsigned long long l, const Integer& n);

		friend giv_all_inlined  Integer& operator *= (Integer& n, const long long l);
		friend giv_all_inlined  Integer& operator *= (Integer& n, const unsigned long long l);
#endif

		//@}

		/*! @name fused add-multiply
		 * @brief
		 * Groups a multiplication adn an addition/division is a single function.
		 * This is usually faster than doing the two operations separately.
		 */
		//@{
		/*! axpy
		 *  <code>res = ax+y</code>.
		 * @param res Integers as in the forumla
		 * @param a Integers as in the forumla
		 * @param x Integers as in the forumla
		 * @param y Integers as in the forumla
		 */
		static giv_all_inlined  Integer& axpy   (Integer& res,
							 const Integer& a, const Integer& x, const Integer& y );
		static giv_all_inlined  Integer& axpy   (Integer& res,
							 const Integer& a, const unsigned long x, const Integer& y );

		/*! axpyin (inplace)
		 *  <code>res += ax</code>.
		 * @param res Integers as in the formula.
		 * @param a Integers as in the formula.
		 * @param x Integers as in the formula.
		 */
		static giv_all_inlined  Integer& axpyin   (Integer& res,
							   const Integer& a, const Integer& x);
		static giv_all_inlined  Integer& axpyin   (Integer& res,
							   const Integer& a, const unsigned long x);

		/*! maxpy
		 *  <code>res = y - ax</code>.
		 * @param res Integers as in the formula.
		 * @param a Integers as in the formula.
		 * @param x Integers as in the formula.
		 * @param y Integers as in the formula.
		 */
		static giv_all_inlined  Integer& maxpy   (Integer& res,
							  const Integer& a, const Integer& x, const Integer& y );
		static giv_all_inlined  Integer& maxpy   (Integer& res,
							  const Integer& a, const unsigned long x, const Integer& y );

		/*! maxpyin
		 *  <code>res -= ax</code>.
		 * @param res Integers as in the formula.
		 * @param a Integers as in the formula.
		 * @param x Integers as in the formula.
		 * @param y Integers as in the formula.
		 */
		static giv_all_inlined  Integer& maxpyin (Integer& res,
							  const Integer& a, const Integer& x );
		static giv_all_inlined  Integer& maxpyin (Integer& res,
							  const Integer& a, const unsigned long x );


		/*! axmy
		 *  <code>res = ax - y</code>.
		 * @param res Integers as in the formula.
		 * @param a Integers as in the formula.
		 * @param x Integers as in the formula.
		 * @param y Integers as in the formula.
		 */
		static giv_all_inlined  Integer& axmy   (Integer& res,
							 const Integer& a, const Integer& x, const Integer& y );
		static giv_all_inlined  Integer& axmy   (Integer& res,
							 const Integer& a, const unsigned long x, const Integer & y );


		/*! axmyin (in place)
		 * <code>res = ax - res</code>.
		 * @param res Integers as in the formula.
		 * @param a Integers as in the formula.
		 * @param x Integers as in the formula.
		 */
		static giv_all_inlined  Integer& axmyin   (Integer& res,
							   const Integer& a, const Integer& x);
		static giv_all_inlined  Integer& axmyin   (Integer& res,
							   const Integer& a, const unsigned long x);
		//@}

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
		//@{
		/*! Division \c q/=d.
		 * @param q quotient
		 * @param d divisor.
		 * @return \c q
		 */
		static giv_all_inlined  Integer& divin (Integer& q, const Integer& d);
		static giv_all_inlined  Integer& divin (Integer& q, const long d);
		static giv_all_inlined  Integer& divin (Integer& q, const unsigned long d);
		/*! Division \c q=n/d.
		 * @param q quotient
		 * @param n dividand.
		 * @param d divisor
		 * @return \c q
		 */
		static giv_all_inlined  Integer& div   (Integer& q, const Integer& n, const Integer& d);
		static giv_all_inlined  Integer& div   (Integer& q, const Integer& n, const long d);
		static giv_all_inlined  Integer& div   (Integer& q, const Integer& n, const int d);
		static giv_all_inlined  Integer& div   (Integer& q, const Integer& n, const unsigned long d);
		/*! Division when \c d divides \c n.
		 * @param q exact quotient
		 * @param n dividend
		 * @param d divisor
		 * @warning if quotient is not exact, the result is not predictable.
		 */
		static giv_all_inlined  Integer& divexact  (Integer& q, const Integer& n, const Integer& d);
		/*! Division when \c d divides \c n.
		 * @param n dividend
		 * @param d divisor
		 * @return  exact quotient \c n/d
		 * @warning if quotient is not exact, the result is not predictable.
		 */
		static giv_all_inlined  Integer  divexact  (const Integer& n, const Integer& d);

		/*! Division operator.
		 * @param d divisor
		 */
		giv_all_inlined Integer  operator /  (const Integer&      d) const;
		giv_all_inlined Integer  operator /  (const unsigned long d) const;
		giv_all_inlined Integer  operator /  (const long          d) const;
		/*! Division operator (inplace).
		 * @param d divisor
		 */
		giv_all_inlined Integer& operator /= (const Integer&      d);
		giv_all_inlined Integer& operator /= (const unsigned long d);
		giv_all_inlined Integer& operator /= (const long          d);
		template<class XXX>
		Integer& operator /=(const XXX& d) {
			return this->operator /= ( (Integer)d );
		}

		// -- operator /
		friend giv_all_inlined  Integer operator / (const int l, const Integer& n);
		friend giv_all_inlined  Integer operator / (const long l, const Integer& n);
		friend giv_all_inlined  Integer operator / (const Integer& n, const int l);
		friend giv_all_inlined  Integer operator / (const Integer& n, const unsigned int l);

		friend giv_all_inlined  Integer& operator /= (Integer& n, const int l);
		friend giv_all_inlined  Integer& operator /= (Integer& n, const long l);
		friend giv_all_inlined  Integer& operator /= (Integer& n, const unsigned int l);

		/*!  Function \c mod (inplace).
		 * \f$ r \gets r \mod n\f$
		 * @param r remainder
		 * @param n modulus
		 */
		static giv_all_inlined  Integer& modin (Integer& r, const Integer& n);
		static giv_all_inlined  Integer& modin (Integer& r, const long n);
		static giv_all_inlined  Integer& modin (Integer& r, const unsigned long n);
		/*!  Function \c mod.
		 * \f$ r \gets n \mod d\f$
		 * @param r remainder
		 * @param n integer
		 * @param d modulus
		 */
		static giv_all_inlined  Integer& mod   (Integer& r, const Integer& n, const Integer& d);
		static giv_all_inlined  Integer& mod   (Integer& r, const Integer& n, const long d);
		static giv_all_inlined  Integer& mod   (Integer& r, const Integer& n, const unsigned long d);

		/*! Euclidean division.
		 * <code> n = d q + r </code>.
		 * Computes both the quotient and the residue (as in quorem).
		 * @param q as in the formula
		 * @param r as in the formula
		 * @param n as in the formula
		 * @param d as in the formula
		 * @return the quotient \c q
		 */
		static giv_all_inlined  Integer& divmod   (Integer& q, Integer& r, const Integer& n, const Integer& d);
		static giv_all_inlined  Integer& divmod   (Integer& q, long& r, const Integer& n, const long d);
		static giv_all_inlined  Integer& divmod   (Integer& q, unsigned long& r, const Integer& n, const unsigned long d);

		/*! @name rounding function
		 * these are the same as the STL ones, except for the signature.
		 * @param res the result
		 * @param n the numerator
		 * @param d the demominator
		 */
		//@{
		static giv_all_inlined  Integer&   ceil (Integer & res, const Integer &n, const Integer & d); // same as std::ceil (n/d)
		static giv_all_inlined  Integer&   floor(Integer & res, const Integer &n, const Integer & d); // same as std::floor(n/d)
		static giv_all_inlined  Integer&   trunc(Integer & res, const Integer &n, const Integer & d); // same as std::trunc(n/d)
		static giv_all_inlined  Integer    ceil (const Integer &n, const Integer & d); // same as std::ceil (n/d)
		static giv_all_inlined  Integer    floor(const Integer &n, const Integer & d); // same as std::floor(n/d)
		static giv_all_inlined  Integer    trunc(const Integer &n, const Integer & d); // same as std::trunc(n/d)
		//@}



		// (FILE gmp++_int_mod.C)
		/*! Modulo operator.
		 * @param n modulus
		 * @return remainder <code> (*this) mod n</code>
		 */
		giv_all_inlined Integer  operator % (const Integer& n) const;
		giv_all_inlined long     operator % (const unsigned long n) const;
		giv_all_inlined long     operator % (const long n) const;
		giv_all_inlined double   operator % (const double n) const;
		short    operator % (const unsigned short n) const
		{
			return (short) ( this->operator % ( (unsigned long)n ) );
		}
		template<class XXX>
		XXX      operator %(const XXX& n) const
		{
			return (XXX)this->operator % ( Integer(n) );
		}

		// -- operator %
		friend giv_all_inlined  Integer operator % (const int l, const Integer& n);
		friend giv_all_inlined  Integer operator % (const long l, const Integer& n);
		friend giv_all_inlined  Integer operator % (const Integer& n, const int l);
		friend giv_all_inlined  Integer operator % (const Integer& n, const unsigned int l);

		friend giv_all_inlined  Integer& operator %= (Integer& n, const int l);
		friend giv_all_inlined  Integer& operator %= (Integer& n, const unsigned int l);


		/*! Modulo operator (inplace).
		 * @param n modulus
		 * @return remainder <code> (*this) <- (*this) mod n</code>
		 */
		giv_all_inlined Integer&  operator %= (const Integer& n);
		giv_all_inlined Integer&  operator %= (const unsigned long n);
		giv_all_inlined Integer&  operator %= (const long n);
#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
		Integer&  operator %= (const long long n) {
			return *this %= (Integer)n;
		}
		Integer&  operator %= (const unsigned long long n) {
			return *this %= (Integer)n;
		}
		giv_all_inlined long long operator % (const long long n) const;
		giv_all_inlined unsigned long long operator % (const unsigned long long n) const;
#endif
		template<class XXX>
		Integer& operator  %=(const XXX& n) {
			return this->operator %= ( (Integer)n );
		}






		// (FILE gmp++_int_gcd.C)
		//------------------------------------- Arithmetic functions
		/*! @name Arithmetic functions */
		//@{
		friend giv_all_inlined  Integer gcd (const Integer& a, const Integer& b);
		friend giv_all_inlined  Integer gcd ( Integer& u, Integer& v,
						      const Integer& a, const Integer& b);
		friend giv_all_inlined  Integer& gcd (Integer& g, const Integer& a, const Integer& b);
		friend giv_all_inlined  Integer& gcd (Integer& g, Integer& u, Integer& v,
						      const Integer& a, const Integer& b);
		// modular inverses
		friend giv_all_inlined  Integer& inv (Integer& u, const Integer& a, const Integer& b);
		friend giv_all_inlined  Integer& invin (Integer& u, const Integer& b);

		friend giv_all_inlined  Integer pp( const Integer& P, const Integer& Q );

		friend giv_all_inlined  Integer& lcm (Integer& g, const Integer& a, const Integer& b);
		friend giv_all_inlined  Integer lcm (const Integer& a, const Integer& b);

		// (FILE gmp++_int_pow.C)
		// - return n^l
		friend giv_all_inlined  Integer& pow(Integer& Res, const Integer& n, const long l);
		friend giv_all_inlined  Integer& pow(Integer& Res, const unsigned long n, const unsigned long l);
		friend giv_all_inlined  Integer& pow(Integer& Res, const Integer& n, const unsigned long l);
		friend giv_all_inlined  Integer& pow(Integer& Res, const Integer& n, const int l)
		{
			return pow(Res, n, (long)l );
		}
		friend giv_all_inlined  Integer& pow(Integer& Res, const Integer& n, const unsigned int l)
		{
			return pow(Res, n, (unsigned long)l );
		}
		friend giv_all_inlined  Integer pow(const Integer& n, const long l);
		friend giv_all_inlined  Integer pow(const Integer& n, const unsigned long l);
		friend giv_all_inlined  Integer pow(const Integer& n, const int l)
		{
			return pow(n, (long)l );
		}
		friend giv_all_inlined  Integer pow(const Integer& n, const unsigned int l)
		{
			return pow(n, (unsigned long)l );
		}

		// - return n^e % m
		friend giv_all_inlined  Integer& powmod(Integer& Res, const Integer& n, const unsigned long e, const Integer& m);
		friend giv_all_inlined  Integer& powmod(Integer& Res, const Integer& n, const long e, const Integer& m);
		friend giv_all_inlined  Integer& powmod(Integer& Res, const Integer& n, const unsigned int e, const Integer& m)
		{
			return powmod(Res, n, (unsigned long)e, m);
		}
		friend giv_all_inlined  Integer& powmod(Integer& Res, const Integer& n, const int e, const Integer& m)
		{
			return powmod(Res, n, (long)e, m);
		}
		friend giv_all_inlined  Integer& powmod(Integer& Res, const Integer& n, const Integer& e, const Integer& m);
		friend giv_all_inlined  Integer powmod(const Integer& n, const unsigned long e, const Integer& m);
		friend giv_all_inlined  Integer powmod(const Integer& n, const long e, const Integer& m);
		friend giv_all_inlined  Integer powmod(const Integer& n, const unsigned int e, const Integer& m)
		{
			return powmod(n, (unsigned long)e, m);
		}
		friend giv_all_inlined  Integer powmod(const Integer& n, const int e, const Integer& m)
		{
			return powmod(n, (long)e, m);
		}
		friend giv_all_inlined  Integer powmod(const Integer& n, const Integer& e, const Integer& m);

		// (FILE gmp++_int_misc.C)
		friend giv_all_inlined  Integer fact ( unsigned long l);

		friend giv_all_inlined  Integer sqrt(const Integer& p);
		friend giv_all_inlined  Integer sqrtrem(const Integer& p, Integer& rem);
		friend giv_all_inlined  Integer& sqrt(Integer& r, const Integer& p);
		friend giv_all_inlined  Integer& sqrtrem(Integer& r, const Integer& p, Integer& rem);


		friend giv_all_inlined  bool root(Integer& q, const Integer&, unsigned int n);
		friend giv_all_inlined  long logp(const Integer& a, const Integer& p) ;
		friend giv_all_inlined  double logtwo(const Integer& a) ;
		friend giv_all_inlined  double naturallog(const Integer& a) ;
		//@}

		//-----------------------------------------Miscellaneous
		/*! @name Miscellaneous.
		*/
		//@{
		friend giv_all_inlined  void swap(Integer& , Integer&);

		friend inline int sign   (const Integer& a)
		{
			return a.priv_sign();
		}

		// (FILE gmp++_int_compare.C)
		// compare to 1 and 0
		friend giv_all_inlined  int isOne(const Integer& a);

		friend giv_all_inlined  int isZero(const Integer& a);


		friend giv_all_inlined  int nonZero(const Integer& a);

		friend giv_all_inlined  int isZero(const short int a);

		friend giv_all_inlined  int isZero(const int a);
		friend giv_all_inlined  int isZero(const long a);
		friend giv_all_inlined  int isZero(const unsigned short int a);
		friend giv_all_inlined  int isZero(const unsigned int a);
		friend giv_all_inlined  int isZero(const unsigned long a);
#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
#if 1 /*  use of C++0x long long integer constant */
		friend giv_all_inlined  int isZero(const unsigned long long a);
#endif
		friend giv_all_inlined  int isZero(const long long a);
#endif
		template<class A,class B>
		static giv_all_inlined bool isleq(const A&a,const B&b)
		{ return a<=b ; }



		// (FILE gmp++_int_misc.C)

		friend giv_all_inlined  int isperfectpower  (const Integer& );

		friend giv_all_inlined  Integer abs(const Integer& n);

		friend giv_all_inlined  Integer& prevprime(Integer&, const Integer& p);
		friend giv_all_inlined  Integer& nextprime(Integer&, const Integer& p);
		friend giv_all_inlined  int probab_prime(const Integer& p);
		friend giv_all_inlined  int probab_prime(const Integer& p, int r);
		friend giv_all_inlined  int jacobi(const Integer& u, const Integer& v) ;
		friend giv_all_inlined  int legendre(const Integer& u, const Integer& v) ;

		Integer& operator++()
		{
			return *this+=1UL;
		} // prefix
		Integer operator++(int)
		{ // postfix
			Integer tmp = *this ;
			++*this;
			return tmp;
		}
		Integer& operator--()
		{
			return *this-=1UL;
		} // prefix
		Integer operator--(int)
		{// postfix
			Integer tmp = *this ;
			--*this;
			return tmp;
		}

		// - return the size in byte
		friend giv_all_inlined  unsigned long length (const Integer& a);
		// - return the size in word.
		giv_all_inlined size_t size() const;
		// - return the size in base B (always exact if B is a power of two)
		giv_all_inlined size_t size_in_base(int B) const;
		// - return the size in bit.
		giv_all_inlined size_t bitsize() const;
		// - return the i-th word of the integer. Word 0 is lowest word.
		giv_all_inlined unsigned long operator[](size_t i) const;

		// -- Convert an Integer to a basic C++ type
		// -- Cast operators
		// -- Cast towards unsigned consider only the absolute value
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
		giv_all_inlined operator unsigned long() const ;
		giv_all_inlined operator long() const ;
#ifndef __GIVARO__DONOTUSE_longlong__
		giv_all_inlined operator unsigned long long() const ;
		giv_all_inlined operator long long() const ;
#endif
		giv_all_inlined operator std::string() const ;
		giv_all_inlined operator float() const ;
		giv_all_inlined operator double() const ;
		giv_all_inlined operator vect_t() const ;
		//@}

		// (FILE gmp++_int_rand.inl)
		/* BB:
		 * if the following functions are NOT static inline, one
		 * can use -Wl,-zmuldefs....
		 */
		//--------------------Random Iterators
		/*! @name Random numbers functions
		*/
		//@{
		static inline void seeding(unsigned long int  s);
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
		template<bool U>
		static inline Integer& random_lessthan (Integer& r, const Integer & m);
		static inline Integer& random_lessthan (Integer& r, const Integer & m) ;
		template<bool U>
		static inline Integer& random_lessthan_2exp (Integer& r, const unsigned long & m);

		static inline Integer& random_lessthan_2exp (Integer& r, const unsigned long & m) ;
		template<bool U>
		static inline Integer random_lessthan_2exp (const unsigned long & m);
		static inline Integer random_lessthan_2exp (const unsigned long & m) ;
		template<bool U>
		static inline Integer& random_lessthan (Integer& r, const unsigned long & m) ;

		static inline Integer& random_lessthan (Integer& r, const unsigned long & m) ;
		template<bool U,class T>
		static inline Integer random_lessthan (const T & m);
		template<class T>
		static inline Integer random_lessthan (const T & m) ;


		/*  random = */
		template<bool U>
		static inline Integer& random_exact_2exp (Integer& r, const unsigned long int & m) ;
		static inline Integer& random_exact_2exp (Integer& r, const unsigned long int & m);


		template<bool U>
		static inline Integer& random_exact (Integer& r, const Integer & s) ;
		static inline Integer& random_exact (Integer& r, const Integer & s) ;

		template<bool U>
		static inline Integer& random_exact (Integer& r, const unsigned long int & m)  ;
		static inline Integer& random_exact (Integer& r, const unsigned long int & m) ;
		template<bool U,class T>
		static inline Integer& random_exact (Integer& r, const T & m) ;
		template<class T>
		static inline Integer& random_exact (Integer& r, const T & m) ;
		template<bool U,class T>
		static inline Integer random_exact (const T & s) ;
		template<class T>
		static inline Integer random_exact (const T & s) ;

		/*  random <.< */
		static inline Integer& random_between (Integer& r, const Integer& m, const Integer&M) ;
		static inline Integer random_between (const Integer& m, const Integer &M) ;
		static inline Integer& random_between_2exp (Integer& r, const unsigned long int& m,
							    const unsigned long int &M) ;
		static inline Integer& random_between (Integer& r, const unsigned long int& m,
						       const unsigned long int &M) ;
		static inline Integer random_between_2exp (const unsigned long int & m,
							   const unsigned long int &M) ;
		static inline Integer random_between (const unsigned long int & m,
						      const unsigned long int &M) ;
		template<class R>
		static inline Integer random_between (const R & m, const R & M) ;
		template<class R>
		static inline Integer & random_between (Integer &r, const R & m, const R & M);


		// useful functions :
		template<bool U,class T>
		static inline Integer& random (Integer& r, const T & m) ;
		template<class T>
		static inline Integer& random (Integer& r, const T & m) ;
		template<bool U,class T>
		static inline Integer random(const T & sz) ;
		template<class T>
		static inline Integer random(const T & sz) ;
		template<bool U>
		static inline Integer random();
		static inline Integer random();
		template<bool U,class T>
		static inline Integer nonzerorandom(const T & sz) ;
		template<bool U,class T>
		static inline Integer& nonzerorandom (Integer& r, const T& size) ;
		template<class T>
		static inline Integer nonzerorandom(const T & sz) ;
		template<class T>
		static inline Integer& nonzerorandom (Integer& r, const T& size) ;
		static inline Integer nonzerorandom() ;
		//@}

		//----------------------------------------------I/O
		/*! @name I/O
		*/
		//@{
		friend giv_all_inlined  std::istream& operator >> (std::istream &i, Integer& n);
		friend giv_all_inlined  std::ostream& operator << (std::ostream &o, const Integer& n);
		friend giv_all_inlined  std::ostream& absOutput (std::ostream &o, const Integer& n);

		friend giv_all_inlined  void importWords(Integer&, size_t, int, int, int, size_t, const void*);

		giv_all_inlined std::ostream& print( std::ostream& o ) const;
		//@}

		int sign() const // BB is this usefull ?
		{
			return priv_sign();
		} // but figure out the friend sign()


		mpz_ptr get_mpz()
		{
			return (mpz_ptr)&gmp_rep;
		}
		mpz_srcptr get_mpz_const() const
		{
			return (mpz_srcptr)&gmp_rep;
		}

		int priv_sign() const // BB no longer protected.
		{
			return mpz_sgn( (mpz_srcptr)&gmp_rep );
		}

	protected:

		typedef __mpz_struct Rep;

		Rep gmp_rep;

		//mpz_ptr get_mpz()
		//{ return (mpz_ptr)&gmp_rep; }

		const Rep* get_rep() const
		{
			return &gmp_rep;
		}


	}; //----------------------------------------------- End of Class Integer

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
