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

#include "recint/ruconvert.h" // For ruint_to_mpz()
#include "recint/rconvert.h" // For rint_to_mpz()

#ifndef __GIVARO_GMPplusplus_H
#warning "you should include gmp++/gmp++.h before gmp++/gmp++_int.h (or prepare for the worse)"
#endif

#ifdef __USE_64_bits__
#    define __USE_GMPPLUSPLUS_SIXTYFOUR__
#endif

#ifdef __USE_ISOC99
#    define __USE_GMPPLUSPLUS_SIXTYFOUR__
#endif

#include <givaro/givcaster.h>

namespace Givaro {

    //------------------------------------------------------ Friend Integer
    // forward declaration.
    class Integer;

    // (FILE gmp++_int_gcd.C)

    /*! @brief Modular inverse.
     * @param a
     * @param b
     * @param[out] u is set to \f$a^{-1}\f$ modulo b
     */
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
    Integer& 	pow(Integer& Res, const Integer& n, const int64_t l);
    Integer& 	pow(Integer& Res, const uint64_t n, const uint64_t l);
    Integer& 	pow(Integer& Res, const Integer& n, const uint64_t l);
    Integer& 	pow(Integer& Res, const Integer& n, const int32_t l) ;
    Integer& 	pow(Integer& Res, const Integer& n, const uint32_t l) ;
    Integer 	pow(const Integer& n, const int64_t l);
    Integer 	pow(const Integer& n, const uint64_t l);
    Integer 	pow(const Integer& n, const int32_t l) ;
    Integer 	pow(const Integer& n, const uint32_t l);

    Integer& 	powmod(Integer& Res, const Integer& n, const uint64_t e, const Integer& m);
    Integer& 	powmod(Integer& Res, const Integer& n, const int64_t e, const Integer& m);
    Integer& 	powmod(Integer& Res, const Integer& n, const uint32_t e, const Integer& m) ;
    Integer& 	powmod(Integer& Res, const Integer& n, const int32_t e, const Integer& m)  ;
    Integer& 	powmod(Integer& Res, const Integer& n, const Integer& e, const Integer& m);
    Integer 	powmod(const Integer& n, const uint64_t e, const Integer& m);
    Integer 	powmod(const Integer& n, const int64_t e, const Integer& m);
    Integer 	powmod(const Integer& n, const uint32_t e, const Integer& m) ;
    Integer 	powmod(const Integer& n, const int32_t e, const Integer& m) ;
    Integer 	powmod(const Integer& n, const Integer& e, const Integer& m);

    // (FILE inline)
    int32_t 		sign   (const Integer& a);

    // (FILE gmp++_int_compare.C)
    int32_t 		compare(const Integer& a, const Integer& b);

    int32_t 		absCompare(const Integer& a, const Integer& b);
    int32_t 		absCompare(const Integer& a, const double b);
    int32_t 		absCompare(const Integer& a, const float b);
    int32_t 		absCompare(const Integer& a, const uint64_t b);
    int32_t 		absCompare(const Integer& a, const unsigned b);
    int32_t 		absCompare(const Integer& a, const int64_t b);
    int32_t 		absCompare(const Integer& a, const int32_t b);
    template<class T>
    int32_t             absCompare( const T a,  const Integer & b)
    {
        return absCompare(b,a);
    }

    int32_t 		isZero (const Integer& a);
    int32_t 		nonZero (const Integer& a);
    int32_t 		isOne  (const Integer& a);
    int32_t 		isMOne  (const Integer& a);

    // (FILE gmp++_int_misc.C)
    Integer 	fact ( uint64_t l);

    Integer 	sqrt(const Integer& p);
    Integer 	sqrtrem(const Integer& p, Integer& rem);
    Integer& 	sqrt(Integer& r, const Integer& p);
    Integer& 	sqrtrem(Integer& r, const Integer& p, Integer& rem);
    bool 		root(Integer& q, const Integer&, uint32_t n);

    int64_t		logp(const Integer& a, const Integer& p) ;
    double 		logtwo(const Integer& a) ;
    double 		naturallog(const Integer& a) ;

    void 		swap(Integer& , Integer&);

    int32_t 		isperfectpower  (const Integer& );

    Integer 	abs(const Integer& n);

    int32_t 		jacobi(const Integer& u, const Integer& v) ;
    int32_t 		legendre(const Integer& u, const Integer& v) ;

    bool            isOdd(const Integer&);

    uint64_t 	length (const Integer& a);
    // (FILE gmp++_int_io.C)

    std::istream& 	operator >> (std::istream &i, Integer& n);
    std::ostream& 	operator << (std::ostream &o, const Integer& n);
    std::ostream& 	absOutput (std::ostream &o, const Integer& n);

    // The following are protected
    // to deal with primality, use Givaro::IntPrimeDom
    namespace Protected {
        void 		importWords(Integer&, size_t, int, int, int, size_t, const void*);
        Integer& 	prevprime(Integer&, const Integer& p);
        Integer& 	nextprime(Integer&, const Integer& p);
        int32_t 		probab_prime(const Integer& p, int32_t r=_GIVARO_ISPRIMETESTS_);
    }


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
        giv_all_inlined Integer(int32_t n = 0);
        //! @overload Givaro::Integer(int)
        giv_all_inlined Integer(int64_t n);
        //! @overload Givaro::Integer(int)
        giv_all_inlined Integer(unsigned char n);
        //! @overload Givaro::Integer(int)
        giv_all_inlined Integer(uint32_t n);
        //! @overload Givaro::Integer(int)
        giv_all_inlined Integer(uint64_t n);
        //! @overload Givaro::Integer(int)
        giv_all_inlined Integer(double n);
        //! @overload Givaro::Integer(int)
        giv_all_inlined Integer(const char *n);
        //! @overload Givaro::Integer(mpz_class)
        giv_all_inlined Integer(const mpz_class& a) {
            mpz_init_set((mpz_ptr)&gmp_rep, a.get_mpz_t()) ;
        }
        //! @overload Givaro::Integer(RecInt::ruint<K>)
        template<size_t K> Integer(const RecInt::ruint<K>& n)
        {
            RecInt::ruint_to_mpz_t(get_mpz(), n);
        }
        //! @overload Givaro::Integer(RecInt::ruint<K>)
        template<size_t K> Integer(const RecInt::rint<K>& n)
        {
            RecInt::rint_to_mpz_t(get_mpz(), n);
        }


        /*! Copy constructor
         * @param n input to be constructed from
         */
        giv_all_inlined Integer(const Integer& n);

        /*! Creates a new Integer from pointers.
         * @param d array
         * @param sz size
         */
        giv_all_inlined Integer(uint64_t* d, int64_t sz);
        /*! Creates a new Integers for a vector of limbs
         * @param v vector of limbs
         */
        giv_all_inlined Integer( const vect_t &v);

        //! destructor
        giv_all_inlined ~Integer();
        ///@}

        // -- type_string
        static const std::string type_string () {
            return "Integer";
        }

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

        //------------------Equalities and inequalities between integers and basic
        // (FILE gmp++_int_compare.C)
        //! @name Comparisons functions.
        ///@{
        /*! Compares two integers.
         * @param a integer
         * @param b integer
         * @return \c 1 if \f$a > b\f$, \c 0 if \f$a = b\f$ and \p -1 otherwise.
         */
        giv_all_inlined friend int32_t compare(const Integer& a, const Integer& b);

        /** Compare the norm of two integers.
         * @param a integer
         * @param b integer
         * @return \c 1 if \f$|a| > |b|\f$, \c 0 if \f$|a| = |b|\f$ and \p -1 otherwise.
         */
        giv_all_inlined friend int32_t absCompare(const Integer& a, const Integer& b);
        /** @overload Integer::absCompare(Integer, Integer) */
        giv_all_inlined friend int32_t absCompare(const Integer& a, const double b);
        /** @overload Integer::absCompare(Integer, Integer) */
        giv_all_inlined friend int32_t absCompare(const Integer& a, const float b);
        /** @overload Integer::absCompare(Integer, Integer) */
        giv_all_inlined friend int32_t absCompare(const Integer& a, const uint64_t b);
        /** @overload Integer::absCompare(Integer, Integer) */
        giv_all_inlined friend int32_t absCompare(const Integer& a, const unsigned b);
        /** @overload Integer::absCompare(Integer, Integer) */
        giv_all_inlined friend int32_t absCompare(const Integer& a, const int64_t b);
        /** @overload Integer::absCompare(Integer, Integer) */
        giv_all_inlined friend int32_t absCompare(const Integer& a, const int32_t b);
        /** @overload Integer::absCompare(Integer, Integer) */
        template<class T>
        giv_all_inlined friend int32_t absCompare( const T a,  const Integer & b) ;

        //! name compare to 1 and 0
        //! @param a
        friend giv_all_inlined  int32_t isOne(const Integer& a);
        friend giv_all_inlined  int32_t isMOne(const Integer& a);



        //! name compare to 1 and 0
        //! @param a
        friend giv_all_inlined  int32_t nonZero(const Integer& a);

        //! name compare to 1 and 0
        //! @param a
        friend giv_all_inlined  int32_t isZero(const Integer& a);
        //! @overload Givaro::isZero(Integer);
        friend giv_all_inlined  int32_t isZero(const int16_t a);
        //! @overload Givaro::isZero(Integer);
        friend giv_all_inlined  int32_t isZero(const int32_t a);
        //! @overload Givaro::isZero(Integer);
        friend giv_all_inlined  int32_t isZero(const int64_t a);
        //! @overload Givaro::isZero(Integer);
        friend giv_all_inlined  int32_t isZero(const uint16_t a);
        //! @overload Givaro::isZero(Integer);
        friend giv_all_inlined  int32_t isZero(const uint32_t a);
        //! @overload Givaro::isZero(Integer);
        friend giv_all_inlined  int32_t isZero(const uint64_t a);
        //! isleq
        //! @param a,b
        template<class A,class B>
        static giv_all_inlined bool isleq(const A&a,const B&b)
        {
            return a<=b ;
        }

        ///@}


        /** @name Comparison operators.
         * @brief Compare with operators.
         */
        ///@{
        /** greater or equal.
         * @param l integer to be compared to
         */
        giv_all_inlined int32_t operator >= (const Integer & l) const;
        /** @overload Integer::operator>=(Integer) */
        giv_all_inlined int32_t operator >= (const int32_t l) const;
        /** @overload Integer::operator>=(Integer) */
        giv_all_inlined int32_t operator >= (const int64_t l) const;
        /** @overload Integer::operator>=(Integer) */
        giv_all_inlined int32_t operator >= (const uint64_t l) const;
        /** @overload Integer::operator>=(Integer) */
        giv_all_inlined int32_t operator >= (const uint32_t l) const;
        /** @overload Integer::operator>=(Integer) */
        giv_all_inlined int32_t operator >= (const double l) const;
        /** @overload Integer::operator>=(Integer) */
        giv_all_inlined int32_t operator >= (const float l) const;

        //! greater or equal.
        /// @param l,n integers to compare
        giv_all_inlined friend int32_t operator >= (uint32_t l, const Integer& n);
        /** @overload Integer::operator>=(unsigned, Integer) */
        giv_all_inlined friend int32_t operator >= (float l, const Integer& n);
        /** @overload Integer::operator>=(unsigned, Integer) */
        giv_all_inlined friend int32_t operator >= (double l, const Integer& n);
        /** @overload Integer::operator>=(unsigned, Integer) */
        giv_all_inlined friend int32_t operator >= (int32_t l, const Integer& n);
        /** @overload Integer::operator>=(unsigned, Integer) */
        giv_all_inlined friend int32_t operator >= (int64_t l, const Integer& n);
        /** @overload Integer::operator>=(unsigned, Integer) */
        giv_all_inlined friend int32_t operator >= (uint64_t l, const Integer& n);


        //! less or equal
        /// @param l integer to be compared to
        giv_all_inlined int32_t operator <= (const Integer & l) const;
        /** @overload Integer::operator<=(Integer) */
        giv_all_inlined int32_t operator <= (const int32_t l) const;
        /** @overload Integer::operator<=(Integer) */
        giv_all_inlined int32_t operator <= (const int64_t l) const;
        /** @overload Integer::operator<=(Integer) */
        giv_all_inlined int32_t operator <= (const uint64_t l) const;
        /** @overload Integer::operator<=(Integer) */
        giv_all_inlined int32_t operator <= (const uint32_t l) const;
        /** @overload Integer::operator<=(Integer) */
        giv_all_inlined int32_t operator <= (const double l) const;
        /** @overload Integer::operator<=(Integer) */
        giv_all_inlined int32_t operator <= (const float l) const;

        //! less or equal
        /// @param l,n integers to compare
        giv_all_inlined friend int32_t operator <= (uint32_t l, const Integer& n);
        /** @overload Integer::operator>=(unsigned, Integer) */
        giv_all_inlined friend int32_t operator <= (float l, const Integer& n);
        /** @overload Integer::operator>=(unsigned, Integer) */
        giv_all_inlined friend int32_t operator <= (double l, const Integer& n);
        /** @overload Integer::operator>=(unsigned, Integer) */
        giv_all_inlined friend int32_t operator <= (int32_t l, const Integer& n);
        /** @overload Integer::operator>=(unsigned, Integer) */
        giv_all_inlined friend int32_t operator <= (int64_t l, const Integer& n);
        /** @overload Integer::operator>=(unsigned, Integer) */
        giv_all_inlined friend int32_t operator <= (uint64_t l, const Integer& n);

        /*! operator != (not equal)
         * @param l integer
         * @return \c 1 iff l == this
         */
        giv_all_inlined int32_t operator != (const Integer & l) const;
        /** @overload Integer::operator!=(Integer) */
        giv_all_inlined int32_t operator != (const int32_t l) const;
        /** @overload Integer::operator!=(Integer) */
        giv_all_inlined int32_t operator != (const int64_t l) const;
        /** @overload Integer::operator!=(Integer) */
        giv_all_inlined int32_t operator != (const uint64_t l) const;
        /** @overload Integer::operator!=(Integer) */
        giv_all_inlined int32_t operator != (const uint32_t l) const;
        /** @overload Integer::operator!=(Integer) */
        giv_all_inlined int32_t operator != (const double l) const;
        /** @overload Integer::operator!=(Integer) */
        giv_all_inlined int32_t operator != (const float l) const;

        /*! operator != (not equal)
         * @param l,n integer
         * @return \c 1 iff l == n
         */
        giv_all_inlined friend int32_t operator != (uint32_t l, const Integer& n);
        /** @overload Integer::operator!=(unsigned, Integer) */
        giv_all_inlined friend int32_t operator != (float l, const Integer& n);
        /** @overload Integer::operator!=(unsigned, Integer) */
        giv_all_inlined friend int32_t operator != (double l, const Integer& n);
        /** @overload Integer::operator!=(unsigned, Integer) */
        giv_all_inlined friend int32_t operator != (int32_t l, const Integer& n);
        /** @overload Integer::operator!=(unsigned, Integer) */
        giv_all_inlined friend int32_t operator != (int64_t l, const Integer& n);
        /** @overload Integer::operator!=(unsigned, Integer) */
        giv_all_inlined friend int32_t operator != (uint64_t l, const Integer& n);


        //! Equality
        /// @param l integer to be compared to
        giv_all_inlined int32_t operator == (const Integer & l) const;
        /** @overload Integer::operator==(Integer) */
        giv_all_inlined int32_t operator == (const int32_t l) const;
        /** @overload Integer::operator==(Integer) */
        giv_all_inlined int32_t operator == (const int64_t l) const;
        /** @overload Integer::operator==(Integer) */
        giv_all_inlined int32_t operator == (const uint64_t l) const;
        /** @overload Integer::operator==(Integer) */
        giv_all_inlined int32_t operator == (const uint32_t l) const;
        /** @overload Integer::operator==(Integer) */
        giv_all_inlined int32_t operator == (const double l) const;
        /** @overload Integer::operator==(Integer) */
        giv_all_inlined int32_t operator == (const float l) const;

        //! Equality
        /// @param l,n integers to compare
        giv_all_inlined friend int32_t operator == (uint32_t l, const Integer& n);
        /** @overload Integer::operator==(unsigned, Integer) */
        giv_all_inlined friend int32_t operator == (float l, const Integer& n);
        /** @overload Integer::operator==(unsigned, Integer) */
        giv_all_inlined friend int32_t operator == (double l, const Integer& n);
        /** @overload Integer::operator==(unsigned, Integer) */
        giv_all_inlined friend int32_t operator == (int32_t l, const Integer& n);
        /** @overload Integer::operator==(unsigned, Integer) */
        giv_all_inlined friend int32_t operator == (int64_t l, const Integer& n);
        /** @overload Integer::operator==(unsigned, Integer) */
        giv_all_inlined friend int32_t operator == (uint64_t l, const Integer& n);


        //! greater (strict)
        /// @param l integer to be compared to
        giv_all_inlined int32_t operator > (const Integer & l) const;
        /** @overload Integer::operator>(Integer) */
        giv_all_inlined int32_t operator > (const int32_t l) const;
        /** @overload Integer::operator>(Integer) */
        giv_all_inlined int32_t operator > (const int64_t l) const;
        /** @overload Integer::operator>(Integer) */
        giv_all_inlined int32_t operator > (const uint64_t l) const;
        /** @overload Integer::operator>(Integer) */
        giv_all_inlined int32_t operator > (const uint32_t l) const;
        /** @overload Integer::operator>(Integer) */
        giv_all_inlined int32_t operator > (const double l) const;
        /** @overload Integer::operator>(Integer) */
        giv_all_inlined int32_t operator > (const float l) const;

        //! greater (strict)
        /// @param l,n integers to compare
        giv_all_inlined friend int32_t operator > (uint32_t l, const Integer& n);
        /** @overload Integer::operator>(unsigned, Integer) */
        giv_all_inlined friend int32_t operator > (float l, const Integer& n);
        /** @overload Integer::operator>(unsigned, Integer) */
        giv_all_inlined friend int32_t operator > (double l, const Integer& n);
        /** @overload Integer::operator>(unsigned, Integer) */
        giv_all_inlined friend int32_t operator > (int32_t l, const Integer& n);
        /** @overload Integer::operator>(unsigned, Integer) */
        giv_all_inlined friend int32_t operator > (int64_t l, const Integer& n);
        /** @overload Integer::operator>(unsigned, Integer) */
        giv_all_inlined friend int32_t operator > (uint64_t l, const Integer& n);

        //! less (strict)
        /// @param l integer to be compared to
        giv_all_inlined int32_t operator < (const Integer & l) const;
        /** @overload Integer::operator<(Integer) */
        giv_all_inlined int32_t operator < (const int32_t l) const;
        /** @overload Integer::operator<(Integer) */
        giv_all_inlined int32_t operator < (const int64_t l) const;
        /** @overload Integer::operator<(Integer) */
        giv_all_inlined int32_t operator < (const uint64_t l) const;
        /** @overload Integer::operator<(Integer) */
        giv_all_inlined int32_t operator < (const uint32_t l) const;
        /** @overload Integer::operator<(Integer) */
        giv_all_inlined int32_t operator < (const double l) const;
        /** @overload Integer::operator<(Integer) */
        giv_all_inlined int32_t operator < (const float l) const;

        //! less (strict)
        /// @param l,n integers to compare
        giv_all_inlined friend int32_t operator < (uint32_t l, const Integer& n);
        /** @overload Integer::operator<(unsigned, Integer) */
        giv_all_inlined friend int32_t operator < (float l, const Integer& n);
        /** @overload Integer::operator<(unsigned, Integer) */
        giv_all_inlined friend int32_t operator < (double l, const Integer& n);
        /** @overload Integer::operator<(unsigned, Integer) */
        giv_all_inlined friend int32_t operator < (int32_t l, const Integer& n);
        /** @overload Integer::operator<(unsigned, Integer) */
        giv_all_inlined friend int32_t operator < (int64_t l, const Integer& n);
        /** @overload Integer::operator<(unsigned, Integer) */
        giv_all_inlined friend int32_t operator < (uint64_t l, const Integer& n);
        ///@}

        // ---------------- Bit logic
        // (FILE gmp++_int_misc.C)

        /*! @name Bit logic operators.  */
        ///@{
        //! @brief XOR (^)
        //! @param a integer
        giv_all_inlined Integer operator^ (const Integer& a) const;   // XOR
        /** @overload Integer::operator^(Integer) */
        giv_all_inlined Integer operator^ (const uint64_t & a) const;
        /** @overload Integer::operator^(Integer) */
        giv_all_inlined Integer operator^ (const uint32_t& a) const;
        //! @brief XOR inplace (^=)
        //! @param a integer
        giv_all_inlined Integer& operator^= (const Integer&a);   // XOR
        /** @overload Integer::operator^=(Integer) */
        giv_all_inlined Integer& operator^= (const uint64_t & a);   // XOR
        /** @overload Integer::operator^=(Integer) */
        giv_all_inlined Integer& operator^= (const uint32_t&a);   // XOR

        //! OR (|)
        //! @param a integer
        giv_all_inlined Integer operator| (const Integer& a ) const;   // OR
        /** @overload Integer::operator|(Integer) */
        giv_all_inlined Integer operator| (const uint64_t & a) const;
        /** @overload Integer::operator|(Integer) */
        giv_all_inlined Integer operator| (const uint32_t& a) const;
        //! OR inplace (|=)
        //! @param a integer
        giv_all_inlined Integer& operator|= (const Integer& a );   // OR
        /** @overload Integer::operator|=(Integer) */
        giv_all_inlined Integer& operator|= (const uint64_t & a );   // OR
        /** @overload Integer::operator|=(Integer) */
        giv_all_inlined Integer& operator|= (const uint32_t& a );   // OR

        //! AND (&)
        //! @param a integer
        giv_all_inlined Integer operator& (const Integer&a) const;   // AND
        /** @overload Integer::operator&(Integer) */
        giv_all_inlined uint32_t operator& (const uint32_t& a) const;
        /** @overload Integer::operator&(Integer) */
        giv_all_inlined uint64_t operator& (const uint64_t & a) const;

        //! AND inplace (&=)
        //! @param a integer
        giv_all_inlined Integer& operator&= (const Integer&a);   // AND
        /** @overload Integer::operator&=(Integer) */
        giv_all_inlined Integer& operator&= (const uint64_t &a);   // AND
        /** @overload Integer::operator&=(Integer) */
        giv_all_inlined Integer& operator&= (const uint32_t&a);   // AND

        //! complement to 1 (~)
        giv_all_inlined Integer operator ~ () const;   // 1 complement

        //! left shift (<<)
        //! @param l shift
        giv_all_inlined Integer operator<< (int32_t l) const; // lshift
        /** @overload Integer::operator<<(int) */
        giv_all_inlined Integer operator<< (int64_t l) const; // lshift
        /** @overload Integer::operator<<(int) */
        giv_all_inlined Integer operator<< (uint32_t l) const; // lshift
        /** @overload Integer::operator<<(int) */
        giv_all_inlined Integer operator<< (uint64_t l) const; // lshift

        //! left shift inplace (<<=)
        //! @param l shift
        giv_all_inlined Integer& operator<<= (int32_t l) ; // lshift
        /** @overload Integer::operator<<=(int) */
        giv_all_inlined Integer& operator<<= (int64_t l) ; // lshift
        /** @overload Integer::operator<<=(Integer) */
        giv_all_inlined Integer& operator<<= (uint32_t l) ; // lshift
        /** @overload Integer::operator<<=(Integer) */
        giv_all_inlined Integer& operator<<= (uint64_t l) ; // lshift

        //! right shift (>>)
        //! @param l shift
        giv_all_inlined Integer operator>> (int32_t l) const; // rshift
        /** @overload Integer::operator>>(int) */
        giv_all_inlined Integer operator>> (int64_t l) const; // rshift
        /** @overload Integer::operator>>(int) */
        giv_all_inlined Integer operator>> (uint32_t l) const; // rshift
        /** @overload Integer::operator>>(int) */
        giv_all_inlined Integer operator>> (uint64_t l) const; // rshift

        //! right shift inplace (>>=)
        //! @param l shift
        giv_all_inlined Integer& operator>>= (int32_t l) ; // rshift
        /** @overload Integer::operator>>=(int) */
        giv_all_inlined Integer& operator>>= (int64_t l) ; // rshift
        /** @overload Integer::operator>>=(int) */
        giv_all_inlined Integer& operator>>= (uint32_t l) ; // rshift
        /** @overload Integer::operator>>=(int) */
        giv_all_inlined Integer& operator>>= (uint64_t l) ; // rshift
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
        static giv_all_inlined  Integer& addin (Integer& res, const int64_t n);
        /** @overload Integer::addin(Integer,Integer) */
        static giv_all_inlined  Integer& addin (Integer& res, const uint64_t n);
        /** @overload Integer::addin(Integer,Integer) */
        static giv_all_inlined  Integer& addin (Integer& res, const int32_t n) {
            return addin(res, (int64_t)n);
        }
        /** @overload Integer::addin(Integer,Integer) */
        static giv_all_inlined  Integer& addin (Integer& res, const uint32_t n) {
            return addin(res,(uint64_t)n);
        }



        /*!  Addition
         * <code>res=n1+n2</code>.
         * @param res as in the formula
         * @param n1 as in the formula
         * @param n2 as in the formula
         */
        static giv_all_inlined  Integer& add   (Integer& res, const Integer& n1, const Integer& n2);
        /** @overload Integer::add(Integer,Integer,Integer) */
        static giv_all_inlined  Integer& add   (Integer& res, const Integer& n1, const int64_t n2);
        /** @overload Integer::add(Integer,Integer,Integer) */
        static giv_all_inlined  Integer& add   (Integer& res, const Integer& n1, const uint64_t n2);
        static giv_all_inlined  Integer& add   (Integer& res, const Integer& n1, const int32_t n2) { return add(res,n1,(int64_t)n2); }
        /** @overload Integer::add(Integer,Integer,Integer) */
        static giv_all_inlined  Integer& add   (Integer& res, const Integer& n1, const uint32_t n2) { return add(res,n1,(uint64_t)n2); }

        // (FILE gmp++_int_sub.C)
        /*!  Substraction (inplace)
         * <code>res-=n</code>.
         * @param res as in the formula
         * @param n as in the formula
         */
        static giv_all_inlined  Integer& subin (Integer& res, const Integer& n);
        /** @overload Integer::subin(Integer,Integer) */
        static giv_all_inlined  Integer& subin (Integer& res, const int64_t n);
        /** @overload Integer::subin(Integer,Integer) */
        static giv_all_inlined  Integer& subin (Integer& res, const uint64_t n);
        /** @overload Integer::subin(Integer,Integer) */
        static giv_all_inlined  Integer& subin (Integer& res, const int32_t n) {
            return subin(res,(int64_t)n); }
        /** @overload Integer::subin(Integer,Integer) */
        static giv_all_inlined  Integer& subin (Integer& res, const uint32_t n) {
            return subin(res,(uint64_t)n); }

        /*!  Substraction
         * <code>res=n1-n2</code>.
         * @param res as in the formula
         * @param n1 as in the formula
         * @param n2 as in the formula
         */
        static giv_all_inlined  Integer& sub   (Integer& res, const Integer& n1, const Integer& n2);
        /** @overload Integer::sub(Integer,Integer,Integer) */
        static giv_all_inlined  Integer& sub   (Integer& res, const Integer& n1, const int64_t n2);
        /** @overload Integer::sub(Integer,Integer,Integer) */
        static giv_all_inlined  Integer& sub   (Integer& res, const Integer& n1, const uint64_t n2);
        static giv_all_inlined  Integer& sub   (Integer& res, const Integer& n1, const int32_t n2) { return sub(res,n1,(uint64_t)n2); }
        /** @overload Integer::sub(Integer,Integer,Integer) */
        static giv_all_inlined  Integer& sub   (Integer& res, const Integer& n1, const uint32_t n2) { return sub(res,n1,(uint64_t)n2); }

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
        static giv_all_inlined  Integer& mulin (Integer& res, const int64_t n);
        /** @overload Integer::mulin(Integer,Integer) */
        static giv_all_inlined  Integer& mulin (Integer& res, const uint64_t n);
        /** @overload Integer::mulin(Integer,Integer) */
        static giv_all_inlined  Integer& mulin (Integer& res, const int32_t n) { return mulin(res,(int64_t)n); }
        /** @overload Integer::mulin(Integer,Integer) */
        static giv_all_inlined  Integer& mulin (Integer& res, const uint32_t n){ return mulin(res,(uint64_t)n); }

        /*! Multiplication
         * <code>res=n1*n2</code>.
         * @param res as in the formula
         * @param n1 as in the formula
         * @param n2 as in the formula
         */
        static giv_all_inlined  Integer& mul   (Integer& res, const Integer& n1, const Integer& n2);
        /** @overload Integer::mul(Integer,Integer,Integer) */
        static giv_all_inlined  Integer& mul   (Integer& res, const Integer& n1, const int64_t n2);
        /** @overload Integer::mul(Integer,Integer,Integer) */
        static giv_all_inlined  Integer& mul   (Integer& res, const Integer& n1, const uint64_t n2);
        /** @overload Integer::mul(Integer,Integer,Integer) */
        static giv_all_inlined  Integer& mul   (Integer& res, const Integer& n1, const int32_t n2) { return mul(res,n1,(int64_t)n2); }
        /** @overload Integer::mul(Integer,Integer,Integer) */
        static giv_all_inlined  Integer& mul   (Integer& res, const Integer& n1, const uint32_t n2) { return mul(res,n1,(uint64_t)n2); }
        ///@}

        //----------------Elementary arithmetic between Integers and basic
        // (FILE gmp++_int_add.C)
        /*! @name Addition, substraction, multiplication operators*/
        ///@{
        /*! operator \c +.
         * @return <code> (*this)+n</code>
         * @param n as in the formula.
         */
        giv_all_inlined Integer  operator + (const Integer& n) const;
        /** @overload Integer::operator+(Integer) */
        giv_all_inlined Integer  operator + (const uint64_t n) const;
        /** @overload Integer::operator+(Integer) */
        giv_all_inlined Integer  operator + (const int64_t n) const;
        /** @overload Integer::operator+(Integer) */
        giv_all_inlined Integer  operator + (const uint32_t n) const { return this->operator+((uint64_t)n); }
        /** @overload Integer::operator+(Integer) */
        giv_all_inlined Integer  operator + (const int32_t n) const { return this->operator+((int64_t)n); }

        /*! operator \c += .
         * @param n asfriend In the formula.
         * @return <code> (*this) += n</code>.
         */
        giv_all_inlined Integer& operator += (const Integer& n);
        /** @overload Integer::operator+=(Integer) */
        giv_all_inlined Integer& operator += (const uint64_t n);
        /** @overload Integer::operator+=(Integer) */
        giv_all_inlined Integer& operator += (const int64_t n);
        /** @overload Integer::operator+=(Integer) */
        giv_all_inlined Integer& operator += (const uint32_t n)  { return this->operator+=((uint64_t)n); }
        /** @overload Integer::operator+=(Integer) */
        giv_all_inlined Integer& operator += (const int32_t n) { return this->operator+=((uint64_t)n); }
        /** @overload Integer::operator+=(Integer) */
        template<class XXX>
        Integer& operator +=(const XXX& n) {
            return this->operator += ( Caster<Integer>(n) );
        }

        // - Friends
        //! operator +.
        //! @param l,n to be added
        friend giv_all_inlined  Integer operator + (const int32_t l, const Integer& n);
        /** @overload friend Integer::operator+(int,Integer) */
        friend giv_all_inlined  Integer operator + (const uint32_t l, const Integer& n);
        /** @overload friend Integer::operator+(int,Integer) */
        friend giv_all_inlined  Integer operator + (const int64_t l, const Integer& n);
        /** @overload friend Integer::operator+(int,Integer) */
        friend giv_all_inlined  Integer operator + (const uint64_t l, const Integer& n);

        // (FILE gmp++_int_sub.C)

        /*! operator \c -.
         * @return <code> (*this)-n</code>
         * @param n as in the formula.
         */
        giv_all_inlined Integer  operator - (const Integer& n) const;
        /** @overload Integer::operator-(Integer) */
        giv_all_inlined Integer  operator - (const uint64_t n) const;
        /** @overload Integer::operator-(Integer) */
        giv_all_inlined Integer  operator - (const int64_t n) const;
        /** @overload Integer::operator+(Integer) */
        giv_all_inlined Integer  operator - (const uint32_t n) const { return this->operator-((uint64_t)n); }
        /** @overload Integer::operator+(Integer) */
        giv_all_inlined Integer  operator - (const int32_t n) const { return this->operator-((int64_t)n); }


        /*! operator \c -= .
         * @param n as in the formula.
         * @return <code> (*this) -= n</code>.
         */
        giv_all_inlined Integer& operator -= (const Integer& n);
        /** @overload Integer::operator-=(Integer) */
        giv_all_inlined Integer& operator -= (const uint64_t n);
        /** @overload Integer::operator-=(Integer) */
        giv_all_inlined Integer& operator -= (const int64_t n);
        /** @overload Integer::operator-=(Integer) */
        giv_all_inlined Integer& operator -= (const uint32_t n) { return this->operator-=((uint64_t)n); }
        /** @overload Integer::operator-=(Integer) */
        giv_all_inlined Integer& operator -= (const int32_t n) { return this->operator-=((int64_t)n); }
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
        friend giv_all_inlined  Integer operator - (const int32_t l, const Integer& n);
        /** @overload friend Integer::operator-(int,Integer) */
        friend giv_all_inlined  Integer operator - (const uint32_t l, const Integer& n);
        /** @overload friend Integer::operator-(int,Integer) */
        friend giv_all_inlined  Integer operator - (const int64_t l, const Integer& n);
        /** @overload friend Integer::operator-(int,Integer) */
        friend giv_all_inlined  Integer operator - (const uint64_t l, const Integer& n);

        // (FILE gmp++_int_mul.C)

        /*! operator \c *.
         * @return <code> (*this)*n</code>
         * @param n as in the formula.
         */
        giv_all_inlined Integer  operator * (const Integer& n) const;
        /** @overload Integer::operator*(Integer) */
        giv_all_inlined Integer  operator * (const uint64_t n) const;
        /** @overload Integer::operator*(Integer) */
        giv_all_inlined Integer  operator * (const int64_t n) const;
        /** @overload Integer::operator*(Integer) */
        giv_all_inlined Integer  operator * (const uint32_t n) const { return this->operator*((uint64_t)n); }
        /** @overload Integer::operator*(Integer) */
        giv_all_inlined Integer  operator * (const int32_t n) const { return this->operator*((int64_t)n); }

        /*! operator \c *= .
         * @param n as in the formula.
         * @return <code> (*this) *= n</code>.
         */
        giv_all_inlined Integer& operator *= (const Integer& n);
        /** @overload Integer::operator*=(Integer) */
        giv_all_inlined Integer& operator *= (const uint64_t n);
        /** @overload Integer::operator*=(Integer) */
        giv_all_inlined Integer& operator *= (const int64_t n);
        /** @overload Integer::operator*=(Integer) */
        giv_all_inlined Integer& operator *= (const uint32_t n) { return this->operator*=((uint64_t)n); }
        /** @overload Integer::operator*=(Integer) */
        giv_all_inlined Integer& operator *= (const int32_t n) { return this->operator*=((int64_t)n); }
        /** @overload Integer::operator*=(Integer) */
        template<class XXX>
        Integer& operator *=(const XXX& n) {
            return this->operator *= ( (Integer)n );
        }

        //! operator *
        //! @param l,n to be multpct
        friend giv_all_inlined  Integer operator * (const int32_t l, const Integer& n);
        /** @overload fried Integer::operator*(int,Integer) */
        friend giv_all_inlined  Integer operator * (const uint32_t l, const Integer& n);
        /** @overload fried Integer::operator*(int,Integer) */
        friend giv_all_inlined  Integer operator * (const int64_t l, const Integer& n);
        /** @overload fried Integer::operator*(int,Integer) */
        friend giv_all_inlined  Integer operator * (const uint64_t l, const Integer& n);

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
                                                 const uint64_t x,
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
                                                   const uint64_t x);

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
                                                  const uint64_t x,
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
                                                  const uint64_t x );

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
                                                 const uint64_t x,
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
                                                   const Integer& a, const uint64_t x);
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
        static giv_all_inlined  Integer& divin (Integer& q, const int64_t d);
        //! @overload Integer::divin(Integer,Integer)
        static giv_all_inlined  Integer& divin (Integer& q, const uint64_t d);

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
                                                const Integer& n, const int64_t d);
        //! @overload Integer::div(Integer,Integer,Integer)
        static giv_all_inlined  Integer& div   (Integer& q,
                                                const Integer& n, const int32_t d);
        //! @overload Integer::div(Integer,Integer,Integer)
        static giv_all_inlined  Integer& div   (Integer& q,
                                                const Integer& n, const uint64_t d);

        /*! Division when \c d divides \c n.
         * @param q exact quotient
         * @param n dividend
         * @param d divisor
         * @warning if quotient is not exact, the result is not predictable.
         */
        static giv_all_inlined  Integer& divexact (Integer& q,
                                                   const Integer& n, const Integer& d);
        static giv_all_inlined  Integer& divexact (Integer& q,
                                                   const Integer& n, const uint64_t & d);
        static giv_all_inlined  Integer& divexact (Integer& q,
                                                   const Integer& n, const int64_t & d);

        /*! Division when \c d divides \c n.
         * @param n dividend
         * @param d divisor
         * @return  exact quotient \c n/d
         * @warning if quotient is not exact, the result is not predictable.
         */
        static giv_all_inlined  Integer  divexact (const Integer& n, const Integer& d);
        static giv_all_inlined  Integer  divexact (const Integer& n, const uint64_t & d);
        static giv_all_inlined  Integer  divexact (const Integer& n, const int64_t & d);

        //! Stuff
        static giv_all_inlined Integer& trem(Integer& r, const Integer &n , const Integer & d);
        static giv_all_inlined Integer& crem(Integer& r, const Integer &n , const Integer & d);
        static giv_all_inlined Integer& frem(Integer& r, const Integer &n , const Integer & d);

        //! Stuff
        static giv_all_inlined Integer& trem(Integer& r, const Integer &n , const uint64_t & d);
        static giv_all_inlined Integer& crem(Integer& r, const Integer &n , const uint64_t & d);
        static giv_all_inlined Integer& frem(Integer& r, const Integer &n , const uint64_t & d);

        //! Stuff
        static giv_all_inlined uint64_t trem(const Integer &n , const uint64_t & d);
        static giv_all_inlined uint64_t crem(const Integer &n , const uint64_t & d);
        static giv_all_inlined uint64_t frem(const Integer &n , const uint64_t & d);

        /*! Division operator.
         * @param d divisor
         */
        giv_all_inlined Integer  operator /  (const Integer&      d) const;
        //! @overload Integer::operator/(Integer)
        giv_all_inlined Integer  operator /  (const uint64_t d) const;
        //! @overload Integer::operator/(Integer)
        giv_all_inlined Integer  operator /  (const int64_t d) const;
        //! @overload Integer::operator/(Integer)
        giv_all_inlined Integer  operator /  (const uint32_t d) const { return this->operator/((uint64_t)d); }
        //! @overload Integer::operator/(Integer)
        giv_all_inlined Integer  operator /  (const int32_t d) const { return this->operator/((int64_t)d); }

        /*! Division operator (inplace).
         * @param d divisor
         */
        giv_all_inlined Integer& operator /= (const Integer&      d);
        //! @overload Integer::operator/=(Integer)
        giv_all_inlined Integer& operator /= (const uint64_t d);
        //! @overload Integer::operator/=(Integer)
        giv_all_inlined Integer& operator /= (const int64_t d);
        //! @overload Integer::operator/=(Integer)
        giv_all_inlined Integer& operator /= (const uint32_t d) { return this->operator/=((uint64_t)d); }
        //! @overload Integer::operator/=(Integer)
        giv_all_inlined Integer& operator /= (const int32_t d) { return this->operator/=((int64_t)d); }
        //! @overload Integer::operator/=(Integer)
        template<class XXX>
        Integer& operator /=(const XXX& d) {
            return this->operator /= ( (Integer)d );
        }

        //! operator /
        friend giv_all_inlined  Integer operator / (const int32_t l, const Integer& n);
        //! @overload Integer::operator/(int,Integer)
        friend giv_all_inlined  Integer operator / (const int64_t l,
                                                    const Integer& n);
        //! operator /
        friend giv_all_inlined  Integer operator / (const uint32_t l, const Integer& n);
        //! @overload Integer::operator/(int,Integer)
        friend giv_all_inlined  Integer operator / (const uint64_t l,
                                                    const Integer& n);

        /*!  Function \c mod (inplace).
         * \f$ r \gets r \mod n\f$
         * @param r remainder
         * @param n modulus
         */
        static giv_all_inlined  Integer& modin (Integer& r, const Integer& n);
        //! @overload Integer::modin(Integer,Integer)
        static giv_all_inlined  Integer& modin (Integer& r, const int64_t n);
        //! @overload Integer::modin(Integer,Integer)
        static giv_all_inlined  Integer& modin (Integer& r, const uint64_t n);
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
                                                const Integer& n, const int64_t d);
        //! @overload Integer::mod(Integer,Integer,Integer)
        static giv_all_inlined  Integer& mod   (Integer& r,
                                                const Integer& n, const uint64_t d);
        //! @overload Integer::mod(Integer,Integer,Integer)
        static giv_all_inlined  Integer& mod   (Integer& r,
                                                const Integer& n, const int32_t d)
        {
            return Integer::mod(r,n,(int64_t)d);
        }

        //! @overload Integer::mod(Integer,Integer,Integer)
        static giv_all_inlined  Integer& mod   (Integer& r,
                                                const Integer& n, const uint32_t d)
        {
            return Integer::mod(r,n,(uint64_t)d);
        }

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
                                                 int64_t & r,
                                                 const Integer& n,
                                                 const int64_t d);
        //! @overload Integer::divmod(Integer,Integer,Integer,Integer)
        static giv_all_inlined  Integer& divmod (Integer& q,
                                                 uint64_t & r,
                                                 const Integer& n,
                                                 const uint64_t d);

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
        giv_all_inlined int64_t     operator % (const uint64_t n) const;
        //! @overload Integer::operator%(Integer);
        giv_all_inlined int64_t     operator % (const int64_t n) const;
        //! @overload Integer::operator%(Integer);
        giv_all_inlined int32_t     operator % (const uint32_t n) const { return (int32_t)this->operator%((uint64_t)n); }
        //! @overload Integer::operator%(Integer);
        giv_all_inlined int32_t     operator % (const int32_t n) const { return (int32_t)this->operator%((int64_t)n); }

        //! @overload Integer::operator%(Integer);
        giv_all_inlined double   operator % (const double n) const;
        //! @overload Integer::operator%(Integer);
        int16_t    operator % (const uint16_t n) const
        {
            return (int16_t) ( this->operator % ( (uint64_t)n ) );
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
        friend giv_all_inlined  Integer operator % (const int64_t l, const Integer& n);
        //! @overload Integer::operator%(int,Integer);
        friend giv_all_inlined  Integer operator % (const uint64_t l, const Integer& n);
        //! @overload Integer::operator%(int,Integer);
        friend giv_all_inlined  Integer operator % (const int32_t l, const Integer& n);

        //! @overload Integer::operator%(int,Integer);
        friend giv_all_inlined  Integer operator % (const uint32_t l, const Integer& n);
        //! @overload Integer::operator%(int,Integer);


        /*! Modulo operator (inplace).
         * @param n modulus
         * @return remainder <code> (*this) <- (*this) mod n</code>
         */
        giv_all_inlined Integer&  operator %= (const Integer& n);
        //! @overload Integer::operator%=(Integer);
        giv_all_inlined Integer&  operator %= (const uint64_t n);
        //! @overload Integer::operator%=(Integer);
        giv_all_inlined Integer&  operator %= (const int64_t n);
        //! @overload Integer::operator%=(Integer);
        giv_all_inlined Integer&  operator %= (const uint32_t n) { return this->operator%=((uint64_t)n); }
        //! @overload Integer::operator%=(Integer);
        giv_all_inlined Integer&  operator %= (const int32_t n) { return this->operator%=((int64_t)n); }
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
        friend giv_all_inlined  Integer& pow(Integer& Res, const Integer& n, const int64_t l);
        //! @overload Integer::pow(Integer,Integer,int64_t)
        friend giv_all_inlined  Integer& pow(Integer& Res, const uint64_t n, const uint64_t l);
        //! @overload Integer::pow(Integer,Integer,int64_t)
        friend giv_all_inlined  Integer& pow(Integer& Res, const Integer& n, const uint64_t l);
        //! @overload Integer::pow(Integer,Integer,int64_t)
        friend giv_all_inlined  Integer& pow(Integer& Res, const Integer& n, const int32_t l)
        {
            return pow(Res, n, (int64_t)l );
        }
        //! @overload Integer::pow(Integer,Integer,int64_t)
        friend giv_all_inlined  Integer& pow(Integer& Res, const Integer& n, const uint32_t l)
        {
            return pow(Res, n, (uint64_t)l );
        }

        //! pow. return \f$n^l\f$
        //!@param n,l
        friend giv_all_inlined  Integer pow(const Integer& n, const int64_t l);
        //! @overload Integer::pow(Integer,int64_t)
        friend giv_all_inlined  Integer pow(const Integer& n, const uint64_t l);
        //! @overload Integer::pow(Integer,int64_t)
        friend giv_all_inlined  Integer pow(const Integer& n, const int32_t l)
        {
            return pow(n, (int64_t)l );
        }
        //! @overload Integer::pow(Integer,int64_t)
        friend giv_all_inlined  Integer pow(const Integer& n, const uint32_t l)
        {
            return pow(n, (uint64_t)l );
        }
        ///@}

        //! modular pow. return \f$n^e \mod  m\f$.
        friend giv_all_inlined  Integer& powmod(Integer& Res, const Integer& n, const uint64_t e, const Integer& m);
        //! @overload Integer::powmod(Integer,Integer,uint64_t,Integer)
        friend giv_all_inlined  Integer& powmod(Integer& Res, const Integer& n, const int64_t e, const Integer& m);
        //! @overload Integer::powmod(Integer,Integer,uint64_t,Integer)
        friend giv_all_inlined  Integer& powmod(Integer& Res, const Integer& n, const uint32_t e, const Integer& m)
        {
            return powmod(Res, n, (uint64_t)e, m);
        }
        //! @overload Integer::powmod(Integer,Integer,uint64_t,Integer)
        friend giv_all_inlined  Integer& powmod(Integer& Res, const Integer& n, const int32_t e, const Integer& m)
        {
            return powmod(Res, n, (int64_t)e, m);
        }
        //! @overload Integer::powmod(Integer,Integer,uint64_t,Integer)
        friend giv_all_inlined  Integer& powmod(Integer& Res, const Integer& n, const Integer& e, const Integer& m);


        //! modular pow. return \f$n^e \mod  m\f$.
        //! @param n,e,m
        friend giv_all_inlined  Integer powmod(const Integer& n, const uint64_t e, const Integer& m);
        //! @overload Integer::powmod(Integer,uint64_t,Integer)
        friend giv_all_inlined  Integer powmod(const Integer& n, const int64_t e, const Integer& m);
        //! @overload Integer::powmod(Integer,uint64_t,Integer)
        friend giv_all_inlined  Integer powmod(const Integer& n, const uint32_t e, const Integer& m)
        {
            return powmod(n, (uint64_t)e, m);
        }
        //! @overload Integer::powmod(Integer,uint64_t,Integer)
        friend giv_all_inlined  Integer powmod(const Integer& n, const int32_t e, const Integer& m)
        {
            return powmod(n, (int64_t)e, m);
        }
        //! @overload Integer::powmod(Integer,uint64_t,Integer)
        friend giv_all_inlined  Integer powmod(const Integer& n, const Integer& e, const Integer& m);

        // (FILE gmp++_int_misc.C)
        //! fact
        //! @param l
        friend giv_all_inlined  Integer fact ( uint64_t l);

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
                                          const Integer&a, uint32_t n);

        //! logs
        //! @param a,p
        friend giv_all_inlined  int64_t logp(const Integer& a, const Integer& p) ;
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
friend inline int32_t sign   (const Integer& a)
{
    return a.priv_sign();
}

//! sign
int32_t sign() const // BB is this usefull ?
{
    return priv_sign();
} // but figure out the friend sign()

//! private sign
int32_t priv_sign() const // BB not protected anymore.
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
mpz_srcptr get_mpz() const
{
    return (mpz_srcptr)&gmp_rep;
}

mpz_srcptr get_mpz_const() const
{
    return (mpz_srcptr)&gmp_rep;
}

/*! returns the number of bytes used to store *this          */
//! @param a
friend giv_all_inlined  uint64_t length (const Integer& a);

/*! returns the number of machine words used to store *this          */
giv_all_inlined size_t size() const;

/*! returns <code>ceil(log_BASE(*this))</code>.        */
giv_all_inlined size_t size_in_base(int32_t B) const;

/*! returns <code>ceil(log_2(*this)) </code>.        */
giv_all_inlined size_t bitsize() const;

/*! return the i-th word of the integer. Word 0 is lowest word.*/
giv_all_inlined uint64_t operator[](size_t i) const;



// (FILE gmp++_int_misc.C)

//! perfect power
friend giv_all_inlined  int32_t isperfectpower  (const Integer& n );

//! absolute value
friend giv_all_inlined  Integer abs(const Integer& n);

//! parity of an integer
friend giv_all_inlined  bool isOdd(const Integer&);

//! @name primes
///@{
friend giv_all_inlined  Integer& Protected::prevprime(Integer&, const Integer& p);
friend giv_all_inlined  Integer& Protected::nextprime(Integer&, const Integer& p);
friend giv_all_inlined  int32_t Protected::probab_prime(const Integer& p, int32_t r);
friend giv_all_inlined  int32_t kronecker(const Integer& u, const Integer& v) ;
friend giv_all_inlined  int32_t jacobi(const Integer& u, const Integer& v) ;
friend giv_all_inlined  int32_t legendre(const Integer& u, const Integer& v) ;
///@}

//! @name Increment/Decrement operators
///@{
Integer& operator++()
{ // prefix
    return *this+=1U;
}
Integer operator++(int)
{ // postfix
    Integer tmp = *this ;
    ++*this;
    return tmp;
}
Integer& operator--()
{// prefix
    return *this-=1U;
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
    return *this!=0U;
}
operator int16_t() const
{
    return (int16_t)(int) *this;
}
operator uint16_t() const
{
    return (uint16_t) (uint32_t) *this;
}
operator unsigned char() const
{
    return (unsigned char) (uint32_t) *this;
}
giv_all_inlined operator uint32_t() const ;
giv_all_inlined operator int32_t() const ;
operator signed char() const
{
    return (signed char) (int) *this;
}
giv_all_inlined operator uint64_t() const ;
giv_all_inlined operator int64_t() const ;
giv_all_inlined operator std::string() const ;
giv_all_inlined operator float() const ;
giv_all_inlined operator double() const ;
giv_all_inlined operator vect_t() const ;
template<size_t K> operator RecInt::ruint<K>() const
{
    RecInt::ruint<K> r;
    return RecInt::mpz_t_to_ruint(r, get_mpz_const());
}
template<size_t K> operator RecInt::rint<K>() const
{
    RecInt::rint<K> r;
    return RecInt::mpz_t_to_rint(r, get_mpz_const());
}
///@}

// (FILE gmp++_int_other.C)
//! Other stuff gmp has (temporary)
///@{
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
static inline void seeding(uint64_t  s);
static inline void seeding(const Integer& s);
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
static inline Integer& random_lessthan_2exp (Integer& r, const uint64_t & m);

static inline Integer& random_lessthan_2exp (Integer& r, const uint64_t & m) ;
template<bool ALWAYSPOSITIVE>
static inline Integer random_lessthan_2exp (const uint64_t & m);
static inline Integer random_lessthan_2exp (const uint64_t & m) ;
template<bool ALWAYSPOSITIVE>
static inline Integer& random_lessthan (Integer& r, const uint64_t & m) ;

static inline Integer& random_lessthan (Integer& r, const uint64_t & m) ;
template<bool ALWAYSPOSITIVE,class T>
static inline Integer random_lessthan (const T & m);
template<class T>
static inline Integer random_lessthan (const T & m) ;


/*  random = */
template<bool ALWAYSPOSITIVE>
static inline Integer& random_exact_2exp (Integer& r, const uint64_t & m) ;
static inline Integer& random_exact_2exp (Integer& r, const uint64_t & m);


template<bool ALWAYSPOSITIVE>
static inline Integer& random_exact (Integer& r, const Integer & s) ;
static inline Integer& random_exact (Integer& r, const Integer & s) ;

template<bool ALWAYSPOSITIVE>
static inline Integer& random_exact (Integer& r, const uint64_t & m)  ;
static inline Integer& random_exact (Integer& r, const uint64_t & m) ;
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
static inline Integer& random_between_2exp (Integer& r, const uint64_t& m,
                                            const uint64_t &M) ;
static inline Integer& random_between (Integer& r, const uint64_t& m,
                                       const uint64_t &M) ;
static inline Integer random_between_2exp (const uint64_t & m,
                                           const uint64_t &M) ;
static inline Integer random_between (const uint64_t & m,
                                      const uint64_t &M) ;
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
friend giv_all_inlined  void Protected::importWords(Integer& x , size_t count, int32_t order,
                                                    int32_t size, int32_t endian, size_t nails,
                                                    const void* op);
/** print32_t integer.
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
static inline int32_t sign (const T a)
{
    // std::cout << ((a>0)-(a<0)) << std::endl;
    return (a>0)-(a<0);
    // if (a == 0) return 0 ;
    // return (a<0)?-1:1 ;
}

template <size_t K>
Integer& Caster(Integer& t, const RecInt::ruint<K>& n) {
    RecInt::ruint_to_mpz_t(t.get_mpz(), n);
    return t;
}

template <size_t K>
Integer& Caster(Integer& t, const RecInt::rint<K>& n) {
    RecInt::rint_to_mpz_t(t.get_mpz(), n);
    return t;
}

template <size_t K>
RecInt::ruint<K>& Caster(RecInt::ruint<K>& t, const Integer&  n) {
    return RecInt::mpz_t_to_ruint(t, n.get_mpz_const());
}

template <size_t K>
RecInt::rint<K>& Caster(RecInt::rint<K>& t, const Integer& n) {
    return RecInt::mpz_t_to_rint(t, n.get_mpz_const());
}


} // Givaro

namespace std {
    // Hashing function for Integer, so that std::set<Integer> can be compiled
    template<>
    struct hash<Givaro::Integer> {
    public:
        size_t operator()(const Givaro::Integer& n) const
        {
            size_t hashValue = 0u;
            auto mpz = n.get_mpz();
            auto size = std::abs(mpz->_mp_size);
            for (auto i = 0; i < size; ++i) {
                hashValue ^= mpz->_mp_d[i];
            }
            return hashValue;
        }
    };
}

// only template code is inlined
#ifdef __GIVARO_INLINE_ALL
#include "gmp++/gmp++_int.C"
#endif
#include "gmp++/gmp++_int_rand.inl"

#endif // __GIVARO_GMPplusplus_integer_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
