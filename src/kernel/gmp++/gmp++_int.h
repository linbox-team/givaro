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


#ifdef __USE_64_bits__
#    define __USE_GMPPLUSPLUS_SIXTYFOUR__
#endif

#ifdef __USE_ISOC99
#    define __USE_GMPPLUSPLUS_SIXTYFOUR__
#endif


//------------------------------------------------------ Friend Integer
// forward declaration.
class Integer;

int 		compare(const Integer& a, const Integer& b);
int 		absCompare(const Integer& a, const Integer& b);
Integer& 	inv (Integer& u, const Integer& a, const Integer& b);
Integer& 	invin (Integer& u, const Integer& b);
Integer 	gcd (const Integer& a, const Integer& b);
Integer 	gcd (const Integer& a, const Integer& b, Integer& u, Integer& v);
Integer& 	gcd (Integer& g, const Integer& a, const Integer& b);
Integer& 	gcd (Integer& g, const Integer& a, const Integer& b, Integer& u, Integer& v);
Integer 	pp( const Integer& P, const Integer& Q );
Integer& 	lcm (Integer& g, const Integer& a, const Integer& b);
Integer 	lcm (const Integer& a, const Integer& b);
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
Integer 	fact ( unsigned long l);
Integer 	sqrt(const Integer& p);
Integer 	sqrtrem(const Integer& p, Integer& rem);
Integer& 	sqrt(Integer& r, const Integer& p);
Integer& 	sqrtrem(Integer& r, const Integer& p, Integer& rem);
bool 		root(Integer& q, const Integer&, unsigned int n);
long 		logp(const Integer& a, const Integer& p) ;
double 		logtwo(const Integer& a) ;
double 		naturallog(const Integer& a) ;
void 		swap(Integer& , Integer&);
int 		sign   (const Integer& a);
int 		isZero (const Integer& a);
int 		isOne  (const Integer& a);
int 		isperfectpower  (const Integer& );
Integer 	abs(const Integer& n);
Integer& 	prevprime(Integer&, const Integer& p);
Integer& 	nextprime(Integer&, const Integer& p);
int 		probab_prime(const Integer& p);
int 		probab_prime(const Integer& p, int r);
int 		jacobi(const Integer& u, const Integer& v) ;
int 		legendre(const Integer& u, const Integer& v) ;
unsigned long 	length (const Integer& a);
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
class Integer
{

public:
	//! vector of limbs (ie a gmp number).
	typedef std::vector<mp_limb_t> vect_t;
	//--------------------------------------cstors & dstors
	/*! @name Constructor/Destructors
	 */
	//@{
	Integer( const std::vector<mp_limb_t>& vect_t );
	Integer(int n = 0);
	Integer(long n);
	Integer(unsigned char n);
	Integer(unsigned int n);
	Integer(unsigned long n);
#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
	Integer(long long n);
	Integer(unsigned long long n);
#endif
	Integer(double d);
	Integer(const char *s);
	Integer(const Integer& n);
	~Integer();
	//@}

	//------------------------------------- predefined null and one
	//! zero
	static const Integer zero;
	//! one
	static const Integer one;

	// -- Assignment and copy operators
	/*! @name Assignment and copy operators
	 */
	//@{
	Integer& operator = (const Integer& n);
	Integer& logcpy(const Integer& n);
	Integer& copy(const Integer& n);
	//@}

	//------------------Equalities and inequalities between integers and longs
	/*! @name (in)equality
	 */
	//@{
	int operator != (const int l) const;
	int operator != (const long l) const;
	int operator != (const unsigned long l) const;

	friend int compare(const Integer& a, const Integer& b);
	friend int absCompare(const Integer& a, const Integer& b);

	int operator > (const int l) const;
	int operator > (const long l) const;
	int operator > (const unsigned long l) const;
	int operator < (const int l) const;
	int operator < (const long l) const;
	int operator < (const unsigned long l) const;
	//@}

	//------------------ Bit logic
	/*!@name Bit logic
	 */
	//@{
	Integer operator^ (const Integer&) const;   // XOR
	Integer operator| (const Integer&) const;   // OR
	Integer operator& (const Integer&) const;   // AND
	unsigned long operator^ (const unsigned long& a) const;
	unsigned long operator| (const unsigned long& a) const;
	unsigned long operator& (const unsigned long& a) const;
	Integer operator ~ () const;   // 1 complement
	Integer& operator^= (const Integer&);   // XOR
	Integer& operator|= (const Integer&);   // OR
	Integer& operator&= (const Integer&);   // AND
	Integer operator<< (int l) const; // lshift
	Integer operator>> (int l) const; // rshift
	Integer operator<< (long l) const; // lshift
	Integer operator>> (long l) const; // rshift
	Integer operator<< (unsigned int l) const; // lshift
	Integer operator>> (unsigned int l) const; // rshift
	Integer operator<< (unsigned long l) const; // lshift
	Integer operator>> (unsigned long l) const; // rshift
	Integer& operator<<= (int l) ; // lshift
	Integer& operator>>= (int l) ; // rshift
	Integer& operator<<= (long l) ; // lshift
	Integer& operator>>= (long l) ; // rshift
	Integer& operator<<= (unsigned int l) ; // lshift
	Integer& operator>>= (unsigned int l) ; // rshift
	Integer& operator<<= (unsigned long l) ; // lshift
	Integer& operator>>= (unsigned long l) ; // rshift
	//@}

	// - Methods
	/*! @name Addition, substraction, multiplication
	 */
	//@{
	/*!  Addition (inplace)
	 * <code>res+=n</code>.
	 * @param res as in the formula
	 * @param n as in the formula
	 */
	static Integer& addin (Integer& res, const Integer& n);
	static Integer& addin (Integer& res, const long n);
	static Integer& addin (Integer& res, const unsigned long n);
	/*!  Addition
	 * <code>res=n1+n2</code>.
	 * @param res as in the formula
	 * @param n1 as in the formula
	 * @param n2 as in the formula
	 */
	static Integer& add   (Integer& res, const Integer& n1, const Integer& n2);
	static Integer& add   (Integer& res, const Integer& n1, const long n2);
	static Integer& add   (Integer& res, const Integer& n1, const unsigned long n2);

	/*!  Substraction (inplace)
	 * <code>res-=n</code>.
	 * @param res as in the formula
	 * @param n as in the formula
	 */
	static Integer& subin (Integer& res, const Integer& n);
	static Integer& subin (Integer& res, const long n);
	static Integer& subin (Integer& res, const unsigned long n);
	/*!  Substraction
	 * <code>res=n1-n2</code>.
	 * @param res as in the formula
	 * @param n1 as in the formula
	 * @param n2 as in the formula
	 */
	static Integer& sub   (Integer& res, const Integer& n1, const Integer& n2);
	static Integer& sub   (Integer& res, const Integer& n1, const long n2);
	static Integer& sub   (Integer& res, const Integer& n1, const unsigned long n2);
	/*!  Negation (inplace)
	 * <code>res=-res</code>.
	 * @param res as in the formula
	 */
	static Integer& negin (Integer& res);
	/*!  Negation
	 * <code>res=-n</code>.
	 * @param n as in the formula
	 * @param res as in the formula
	 */
	static Integer& neg   (Integer& res, const Integer& n);

	/*!  Multiplication (inplace)
	 * <code>res*=n</code>.
	 * @param res as in the formula
	 * @param n as in the formula
	 */
	static Integer& mulin (Integer& res, const Integer& n);
	static Integer& mulin (Integer& res, const long n);
	static Integer& mulin (Integer& res, const unsigned long n);
	/*! Multiplication
	 * <code>res=n1*n2</code>.
	 * @param res as in the formula
	 * @param n1 as in the formula
	 * @param n2 as in the formula
	 */
	static Integer& mul   (Integer& res, const Integer& n1, const Integer& n2);
	static Integer& mul   (Integer& res, const Integer& n1, const long n2);
	static Integer& mul   (Integer& res, const Integer& n1, const unsigned long n2);

	//----------------Elementary arithmetic between Integers & longs
	/*! operator \c +.
	 * @return <code> (*this)+n</code>
	 * @param n as in the formula.
	 */
	Integer  operator + (const Integer& n) const;
	Integer  operator + (const unsigned long n) const;
	Integer  operator + (const long n) const;
	/*! operator \c += .
	 * @param n as in the formula.
	 * @return <code> (*this) += n</code>.
	 */
	Integer& operator += (const Integer& n);
	Integer& operator += (const unsigned long n);
	Integer& operator += (const long n);
	template<class XXX>
	Integer& operator +=(const XXX& n) { return this->operator += ( (Integer)n ); }

	/*! operator \c -.
	 * @return <code> (*this)-n</code>
	 * @param n as in the formula.
	 */
	Integer  operator - (const Integer& n) const;
	Integer  operator - (const unsigned long n) const;
	Integer  operator - (const long n) const;
	/*! operator \c -= .
	 * @param n as in the formula.
	 * @return <code> (*this) -= n</code>.
	 */
	Integer& operator -= (const Integer& n);
	Integer& operator -= (const unsigned long n);
	Integer& operator -= (const long n);
	template<class XXX>
	Integer& operator -=(const XXX& n) { return this->operator -= ( (Integer)n ); }

	/*! Opposite.
	 * \return <code>-(*this)</code>.
	 */
	Integer  operator -() const;

	/*! operator \c *.
	 * @return <code> (*this)*n</code>
	 * @param n as in the formula.
	 */
	Integer  operator * (const Integer& n) const;
	Integer  operator * (const unsigned long n) const;
	Integer  operator * (const long n) const;
	/*! operator \c *= .
	 * @param n as in the formula.
	 * @return <code> (*this) *= n</code>.
	 */
	Integer& operator *= (const Integer& n);
	Integer& operator *= (const unsigned long n);
	Integer& operator *= (const long n);
	template<class XXX>
	Integer& operator *=(const XXX& n) { return this->operator *= ( (Integer)n ); }
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
	static Integer& axpy   (Integer& res,
				const Integer& a, const Integer& x, const Integer& y );
	static Integer& axpy   (Integer& res,
				const Integer& a, const unsigned long x, const Integer& y );

	/*! axpyin (inplace)
	 *  <code>res += ax</code>.
	 * @param res Integers as in the formula.
	 * @param a Integers as in the formula.
	 * @param x Integers as in the formula.
	 */
	static Integer& axpyin   (Integer& res,
				  const Integer& a, const Integer& x);
	static Integer& axpyin   (Integer& res,
				  const Integer& a, const unsigned long x);

	/*! maxpy
	 *  <code>res = y - ax</code>.
	 * @param res Integers as in the formula.
	 * @param a Integers as in the formula.
	 * @param x Integers as in the formula.
	 * @param y Integers as in the formula.
	 */
	static Integer& maxpy   (Integer& res,
				 const Integer& a, const Integer& x, const Integer& y );
	static Integer& maxpy   (Integer& res,
				 const Integer& a, const unsigned long x, const Integer& y );

	/*! maxpyin
	 *  <code>res -= ax</code>.
	 * @param res Integers as in the formula.
	 * @param a Integers as in the formula.
	 * @param x Integers as in the formula.
	 * @param y Integers as in the formula.
	 */
	static Integer& maxpyin (Integer& res,
				 const Integer& a, const Integer& x );
	static Integer& maxpyin (Integer& res,
				 const Integer& a, const unsigned long x );


	/*! axmy
	*  <code>res = ax - y</code>.
	 * @param res Integers as in the formula.
	 * @param a Integers as in the formula.
	 * @param x Integers as in the formula.
	 * @param y Integers as in the formula.
	*/
	static Integer& axmy   (Integer& res,
				const Integer& a, const Integer& x, const Integer& y );
	static Integer& axmy   (Integer& res,
				const Integer& a, const unsigned long x, const Integer & y );


	/*! axmyin (in place)
	 * <code>res = ax - res</code>.
	 * @param res Integers as in the formula.
	 * @param a Integers as in the formula.
	 * @param x Integers as in the formula.
	 */
	static Integer& axmyin   (Integer& res,
				  const Integer& a, const Integer& x);
	static Integer& axmyin   (Integer& res,
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
	//@{
	/*! Division \c q/=d.
	 * @param q quotient
	 * @param d divisor.
	 * @return \c q
	 */
	static Integer& divin (Integer& q, const Integer& d);
	static Integer& divin (Integer& q, const long d);
	static Integer& divin (Integer& q, const unsigned long d);
	/*! Division \c q=n/d.
	 * @param q quotient
	 * @param n dividand.
	 * @param d divisor
	 * @return \c q
	 */
	static Integer& div   (Integer& q, const Integer& n, const Integer& d);
	static Integer& div   (Integer& q, const Integer& n, const long d);
	static Integer& div   (Integer& q, const Integer& n, const unsigned long d);
	/*! Division when \c d divides \c n.
	 * @param q exact quotient
	 * @param n dividend
	 * @param d divisor
	 * @warning if quotient is not exact, the result is not predictable.
	 */
	static Integer& divexact  (Integer& q, const Integer& n, const Integer& d);
	/*! Division when \c d divides \c n.
	 * @param n dividend
	 * @param d divisor
	 * @return  exact quotient \c n/d
	 * @warning if quotient is not exact, the result is not predictable.
	 */
	static Integer  divexact  (const Integer& n, const Integer& d);

	/*! Division operator.
	 * @param d divisor
	 */
	Integer  operator /  (const Integer&      d) const;
	Integer  operator /  (const unsigned long d) const;
	Integer  operator /  (const long          d) const;
	/*! Division operator (inplace).
	 * @param d divisor
	 */
	Integer& operator /= (const Integer&      d);
	Integer& operator /= (const unsigned long d);
	Integer& operator /= (const long          d);
	template<class XXX>
	Integer& operator /=(const XXX& d) { return this->operator /= ( (Integer)d ); }

	/*!  Function \c mod (inplace).
	 * \f$ r \gets r \mod n\f$
	 * @param r remainder
	 * @param n modulus
	 */
	static Integer& modin (Integer& r, const Integer& n);
	static Integer& modin (Integer& r, const long n);
	static Integer& modin (Integer& r, const unsigned long n);
	/*!  Function \c mod.
	 * \f$ r \gets n \mod d\f$
	 * @param r remainder
	 * @param n integer
	 * @param d modulus
	 */
	static Integer& mod   (Integer& r, const Integer& n, const Integer& d);
	static Integer& mod   (Integer& r, const Integer& n, const long d);
	static Integer& mod   (Integer& r, const Integer& n, const unsigned long d);

	/*! Euclidean division.
	 * <code> n = d q + r </code>.
	 * @param q as in the formula
	 * @param r as in the formula
	 * @param n as in the formula
	 * @param d as in the formula
	 * @return the quotient \c q
	 */
	static Integer& divmod   (Integer& q, Integer& r, const Integer& n, const Integer& d);
	static Integer& divmod   (Integer& q, long& r, const Integer& n, const long d);
	static Integer& divmod   (Integer& q, unsigned long& r, const Integer& n, const unsigned long d);

	/*! Modulo operator.
	 * @param n modulus
	 * @return remainder <code> (*this) mod n</code>
	 */
	Integer  operator % (const Integer& n) const;
	long     operator % (const unsigned long n) const;
	long     operator % (const long n) const;
	double   operator % (const double n) const;
	short    operator % (const unsigned short n) const { return (short) ( this->operator % ( (unsigned long)n ) ); }
	template<class XXX>
	XXX      operator %(const XXX& n) const { return (XXX)this->operator % ( Integer(n) ); }

	/*! Modulo operator (inplace).
	 * @param n modulus
	 * @return remainder <code> (*this) <- (*this) mod n</code>
	 */
	Integer&  operator %= (const Integer& n);
	Integer&  operator %= (const unsigned long n);
	Integer&  operator %= (const long n);
#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
	Integer&  operator %= (const long long n) { return *this %= (Integer)n; }
	Integer&  operator %= (const unsigned long long n) { return *this %= (Integer)n; }
	long long operator % (const long long n) const;
	unsigned long long operator % (const unsigned long long n) const;
#endif
	template<class XXX>
	Integer& operator  %=(const XXX& n) { return this->operator %= ( (Integer)n ); }
	//@}


	//------------------------------------- Arithmetic functions
	/*! @name Arithmetic functions */
	//@{
	friend Integer gcd (const Integer& a, const Integer& b);
	friend Integer gcd (const Integer& a, const Integer& b,
			    Integer& u, Integer& v);
	friend Integer& gcd (Integer& g, const Integer& a, const Integer& b);
	friend Integer& gcd (Integer& g, const Integer& a, const Integer& b,
			     Integer& u, Integer& v);
        // modular inverses
	friend Integer& inv (Integer& u, const Integer& a, const Integer& b);
	friend Integer& invin (Integer& u, const Integer& b);

	friend Integer pp( const Integer& P, const Integer& Q );

	friend Integer& lcm (Integer& g, const Integer& a, const Integer& b);
	friend Integer lcm (const Integer& a, const Integer& b);

	// - return n^l
	friend Integer& pow(Integer& Res, const Integer& n, const long l);
	friend Integer& pow(Integer& Res, const unsigned long n, const unsigned long l);
	friend Integer& pow(Integer& Res, const Integer& n, const unsigned long l);
	friend Integer& pow(Integer& Res, const Integer& n, const int l) { return pow(Res, n, (long)l ); }
	friend Integer& pow(Integer& Res, const Integer& n, const unsigned int l) { return pow(Res, n, (unsigned long)l ); }
	friend Integer pow(const Integer& n, const long l);
	friend Integer pow(const Integer& n, const unsigned long l);
	friend Integer pow(const Integer& n, const int l) { return pow(n, (long)l ); }
	friend Integer pow(const Integer& n, const unsigned int l) { return pow(n, (unsigned long)l ); }

	// - return n^e % m
	friend Integer& powmod(Integer& Res, const Integer& n, const unsigned long e, const Integer& m);
	friend Integer& powmod(Integer& Res, const Integer& n, const long e, const Integer& m);
	friend Integer& powmod(Integer& Res, const Integer& n, const unsigned int e, const Integer& m) { return powmod(Res, n, (unsigned long)e, m); }
	friend Integer& powmod(Integer& Res, const Integer& n, const int e, const Integer& m)  { return powmod(Res, n, (long)e, m); }
	friend Integer& powmod(Integer& Res, const Integer& n, const Integer& e, const Integer& m);
	friend Integer powmod(const Integer& n, const unsigned long e, const Integer& m);
	friend Integer powmod(const Integer& n, const long e, const Integer& m);
	friend Integer powmod(const Integer& n, const unsigned int e, const Integer& m) { return powmod(n, (unsigned long)e, m); }
	friend Integer powmod(const Integer& n, const int e, const Integer& m)  { return powmod(n, (long)e, m); }
	friend Integer powmod(const Integer& n, const Integer& e, const Integer& m);

	friend Integer fact ( unsigned long l);

	friend Integer sqrt(const Integer& p);
	friend Integer sqrtrem(const Integer& p, Integer& rem);
	friend Integer& sqrt(Integer& r, const Integer& p);
	friend Integer& sqrtrem(Integer& r, const Integer& p, Integer& rem);


	friend bool root(Integer& q, const Integer&, unsigned int n);
	friend long logp(const Integer& a, const Integer& p) ;
	friend double logtwo(const Integer& a) ;
	friend double naturallog(const Integer& a) ;
	//@}

	//-----------------------------------------Miscellaneous
	/*! @name Miscellaneous.
	 */
	//@{
	friend void swap(Integer& , Integer&);

	friend inline int sign   (const Integer& a);
	friend inline int isZero (const Integer& a);
	friend inline int isOne  (const Integer& a);
	friend int isperfectpower  (const Integer& );

	friend Integer abs(const Integer& n);

	friend Integer& prevprime(Integer&, const Integer& p);
	friend Integer& nextprime(Integer&, const Integer& p);
	friend int probab_prime(const Integer& p);
	friend int probab_prime(const Integer& p, int r);
	friend int jacobi(const Integer& u, const Integer& v) ;
	friend int legendre(const Integer& u, const Integer& v) ;

	Integer& operator++() { return *this+=1UL; } // prefix
	Integer operator++(int)
	{ // postfix
	       Integer tmp = *this ;
	       ++*this;
       	       return tmp;
	}
	Integer& operator--() { return *this-=1UL; } // prefix
	Integer operator--(int)
	{// postfix
	       Integer tmp = *this ;
	       --*this;
       	       return tmp;
	}

	// - return the size in byte
	friend inline unsigned long length (const Integer& a);
	// - return the size in word.
	size_t size() const;
	// - return the size in base B (always exact if B is a power of two)
	size_t size_in_base(int B) const;
	// - return the size in bit.
	size_t bitsize() const;
	// - return the i-th word of the integer. Word 0 is lowest word.
	unsigned long operator[](size_t i) const;

	// -- Convert an Integer to a basic C++ type
	// -- Cast operators
	operator bool() const { return *this!=0UL; }
	operator short() const { return (short)(int) *this; }
	operator unsigned short() const { return (unsigned short) (unsigned int) *this; }
	operator unsigned char() const { return (unsigned char) (unsigned int) *this; }
	operator unsigned int() const ;
	operator int() const ;
	operator signed char() const { return (signed char) (int) *this; }
	operator unsigned long() const ;
	operator long() const ;
#ifndef __GIVARO__DONOTUSE_longlong__
	operator unsigned long long() const ;
	operator long long() const ;
#endif
	operator std::string() const ;
	operator float() const ;
	operator double() const ;
	operator vect_t() const ;
	//@}

	//--------------------Random Iterators
	/*! @name Random numbers functions
	 */
	//@{
	static void seeding(unsigned long int  s);
	static void seeding(Integer s);
	static void seeding();

#ifdef __GMP_PLUSPLUS__
	static gmp_randclass& randstate();
#else
	static __gmp_randstate_struct intializerandstate();
	static __gmp_randstate_struct* randstate();
#endif
	static bool RandBool()  ;
	/*  random <= */
	template<bool U>
	static Integer& random_lessthan (Integer& r, const Integer & m);
	static Integer& random_lessthan (Integer& r, const Integer & m)
	{return random_lessthan<true>(r,m);}
	template<bool U>
	static Integer& random_lessthan_2exp (Integer& r, const unsigned long & m);
	static Integer& random_lessthan_2exp (Integer& r, const unsigned long & m)
	{ return random_lessthan_2exp<true>(r,m);}
	template<bool U>
	static Integer random_lessthan_2exp (const unsigned long & m);
	static Integer random_lessthan_2exp (const unsigned long & m)
	{ return random_lessthan_2exp<true>(m);}
	template<bool U>
	static Integer& random_lessthan (Integer& r, const unsigned long & m) ;
	static Integer& random_lessthan (Integer& r, const unsigned long & m)
	{ return random_lessthan<true>(r,m);}
	template<bool U,class T>
	static Integer random_lessthan (const T & m);
	template<class T>
	static Integer random_lessthan (const T & m)
	{ return random_lessthan<true>(m);}


	/*  random = */
	template<bool U>
	static Integer& random_exact (Integer& r, const Integer & s) ;
	static Integer& random_exact (Integer& r, const Integer & s)
	{ return random_exact<true>(r,s); }
	template<bool U>
	static Integer& random_exact_2exp (Integer& r, const unsigned long int & m) ;
	static Integer& random_exact_2exp (Integer& r, const unsigned long int & m)
	{return random_exact_2exp<true>(r,m);}
	template<bool U>
	static Integer& random_exact (Integer& r, const unsigned long int & m)  ;
	static Integer& random_exact (Integer& r, const unsigned long int & m)
	{return random_exact<true>(r,m);}
	template<bool U,class T>
	static Integer& random_exact (Integer& r, const T & m)
	{ return random_exact<U>(r,static_cast<unsigned long int>(m)); }
	template<class T>
	static Integer& random_exact (Integer& r, const T & m)
	{ return random_exact(r,static_cast<unsigned long int>(m)); }
	template<bool U,class T>
	static Integer random_exact (const T & s) ;
	template<class T>
	static Integer random_exact (const T & s)
	{ return random_exact<true>(s) ; }

	/*  random <.< */
	static Integer& random_between (Integer& r, const Integer& m, const Integer&M) ;
	static Integer random_between (const Integer& m, const Integer &M) ;
	static Integer& random_between_2exp (Integer& r, const unsigned long int& m, const unsigned long int &M) ;
	static Integer& random_between (Integer& r, const unsigned long int& m, const unsigned long int &M) ;
	static Integer random_between_2exp (const unsigned long int & m, const unsigned long int &M) ;
	static Integer random_between (const unsigned long int & m, const unsigned long int &M) ;

	template<class R>
	static Integer random_between (const R & m, const R & M)
	{ return random_between(static_cast<unsigned long int>(m), static_cast<unsigned long int>(M)); }

	template<class R>
	static Integer & random_between (Integer &r, const R & m, const R & M)
	{ return random_between(r,static_cast<unsigned long int>(m), static_cast<unsigned long int>(M)); }


	// useful functions :
	template<bool U,class T>
	static Integer& random (Integer& r, const T & m) ;
	template<class T>
	static Integer& random (Integer& r, const T & m)
	{return random<true>(r,m);}
	template<bool U,class T>
	static Integer random(const T & sz) ;
	template<class T>
	static Integer random(const T & sz)
	{ return random<true>(sz);}
	template<bool U>
	static Integer random();
	static Integer random();
	template<bool U,class T>
	static Integer nonzerorandom(const T & sz) ;
	template<bool U,class T>
	static Integer& nonzerorandom (Integer& r, const T& size) ;
	template<class T>
	static Integer nonzerorandom(const T & sz)
	{ return nonzerorandom<true>(sz); }
	template<class T>
	static Integer& nonzerorandom (Integer& r, const T& size)
	{ return nonzerorandom<true>(r,size); }
	static Integer nonzerorandom()
	{
		Integer rez = Integer::nonzerorandom(sizeof(mp_limb_t)*8) ;
		// if (!U) if (Integer::RandBool()) negin(rez);
		return rez;
	}
	//@}

	//----------------------------------------------I/O
	/*! @name I/O
	 */
	//@{
	friend std::istream& operator >> (std::istream &i, Integer& n);
	friend std::ostream& operator << (std::ostream &o, const Integer& n);
	friend std::ostream& absOutput (std::ostream &o, const Integer& n);

	friend void importWords(Integer&, size_t, int, int, int, size_t, const void*);

	std::ostream& print( std::ostream& o ) const;
	//@}

	int sign() const {return priv_sign(); } // but figure out the friend sign()

	mpz_ptr get_mpz() {return (mpz_ptr)&gmp_rep;}

protected:

	typedef __mpz_struct Rep;

	Rep gmp_rep;

	int priv_sign() const;
	//mpz_ptr get_mpz() {return (mpz_ptr)&gmp_rep;}
	const Rep* get_rep() const { return &gmp_rep; }

	// -- Creates a new Integer from a size sz and a array of unsigned long d
	Integer(unsigned long* d, long size);

}; //----------------------------------------------- End of Class Integer



#include "gmp++/gmp++_int.inl"

#endif // __GIVARO_GMPplusplus_integer_H
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen
