// ==========================================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// file: gfq.h
// Time-stamp: <25 Aug 22 14:28:13 Jean-Guillaume.Dumas@imag.fr>
// date: 1999
// version:
// author: Jean-Guillaume.Dumas

/*! @file gfq.h
 * @ingroup zpz
 * @brief   Arithmetic on GF(p^k), with p a prime number less than 2^16.
 */

#ifndef __GIVARO_gfq1_H
#define __GIVARO_gfq1_H

#include "givaro/givconfig.h"
#include "givaro/givinteger.h"
#include "givaro/givranditer.h"
#include "givaro/givpoly1factor.h"
#include "givaro/modular-implem.h"

#include <string>
#include <vector>

#ifdef __GIVARO_COUNT__
#include <iostream>
#endif

namespace Givaro {

    template<class Size_store> class GFqDom;

        // Simplified alias
    template<class Size_store=int64_t> using GFq=GFqDom<Size_store>;



//! class GFqDom
template<class Size_store> class GFqDom {
protected:
	typedef typename Signed_Trait<Size_store>::unsigned_type UTT;
	typedef typename Signed_Trait<Size_store>::signed_type TT;
		// Forces signed storage
	typedef TT Rep;
	typedef typename std::vector<UTT>::size_type  UT  ;
public:
	Rep zero;
	Rep one;
    Rep mOne;
protected:
	UTT _characteristic;	// Field Characteristic (p)
	UTT _exponent;		// Extension degree (k)
	UTT _irred;		// Irreducible polynomial in p-adic
	UTT _q;			// p^k
	UTT _qm1;		// p^k-1

	// G is a generator of GF(q)
	// p is GF(q)'s characteristic
	// log2pol[ i ] = G^i(p)
	// pol2log[ j ] = i such that log2pol[i] = j
	// plus1[i] = k such that G^i + 1 = G^k
	std::vector<UTT> _log2pol;
	std::vector<UTT> _pol2log;
	std::vector<TT> _plus1;

	// Floating point representations
	double _dcharacteristic;

public:

    UTT zech2padic(UTT x) const { return _log2pol[x]; };
    UTT padic2zech(UTT x) const { return _pol2log[x]; };

public:
	typedef GFqDom<TT> Self_t;
	typedef Rep Element;
	typedef Element* Element_ptr ;
	typedef const Element* ConstElement_ptr;


	//     class Element {
	//     public:
	//         mutable Rep _value;

	//         Element() {}
	//     };

	typedef UTT Residu_t;

	// ----- Representation of vector of the Element
	typedef Rep* Array;
	typedef const Rep* constArray;

	GFqDom(): zero(0), one(1), mOne(-1), _log2pol(0), _pol2log(0),_plus1(0) {}

        // Automatic construction
	GFqDom( const UTT P, const UTT e = 1);

        // Construction with prescribed irreducible polynomial
        //   coefficients of the vector should be integers-like
        //   there will be a call to this->init to build the
        //   representation of the irreducible polynomial
    template<typename Vector>
    GFqDom(const UTT P, const UTT e, const Vector& modPoly);

        // Construction with prescribed irreducible polynomial
        //   and with prescribed generator polynomial
        //   coefficients of the vector should be integers-like
        //   there will be a call to this->init to build the
        //   representation of both polynomials
    template<typename Vector>
    GFqDom( const UTT P, const UTT e, const Vector& modPoly, const Vector& genPoly);

    GFqDom( const GFqDom<TT>& F)
            : zero(F.zero),
              one(F.one),
              mOne(F.mOne),
              _characteristic(F._characteristic),
              _exponent(F._exponent),
              _irred(F._irred),
              _q(F._q),
              _qm1(F._qm1),
              _log2pol(F._log2pol),
              _pol2log(F._pol2log),
              _plus1(F._plus1),
              _dcharacteristic(F._dcharacteristic)
        {}

	// Allows to choose the randomization
	// and therefore the field generator
	//     template<class RandIter >
	//     GFqDom(RandIter& g, const UTT P, const UTT e = 1);

	// Destructor
	//    ~GFqDom() {};


    GFqDom<TT>& operator=( const GFqDom<TT>& F) {
        this->zero = F.zero;
        this->one = F.one;
        this->mOne = F.mOne;
        this->_characteristic = F._characteristic;
        this->_dcharacteristic = F._dcharacteristic;
        this->_exponent = F._exponent;
        this->_irred = F._irred;
        this->_q = F._q;
        this->_qm1 = F._qm1;
        this->_log2pol = F._log2pol;
        this->_pol2log = F._pol2log;
        this->_plus1 = F._plus1;
        return *this;
    }

	// ----- Accessors
	inline Element minElement() const { return zero; }
	inline Element maxElement() const { return _qm1; }

	// Access to the modulus, characteristic, size, exponent
    UTT residu() const;
    UTT characteristic() const;
    Integer& characteristic(Integer& p) const {
        return p=characteristic();
    }
    uint64_t& characteristic(uint64_t& p) const {
        return p=(uint64_t)_characteristic;
    }

        // -- maxCardinality
        // -- mostly the limit is store 3 p-large tables
        // -- with int64_t it means 96 GB of RAM
        // -- and about 960 CPU seconds to populate them
    __GIVARO_CONDITIONAL_TEMPLATE(S = TT, (sizeof(S) == sizeof(int32_t)))
    static inline Residu_t maxCardinality() { return 65536u; }
    __GIVARO_CONDITIONAL_TEMPLATE(S = TT, (sizeof(S) == sizeof(int64_t)))
    static inline Residu_t maxCardinality() { return UINT64_C(4294967296); }


    static inline Residu_t minCardinality() { return 2; }

    UTT cardinality() const;
    template<typename T> T& cardinality(T& p) const { return p = T(_q); }
    UTT size() const;
    UTT exponent() const;
	// Internal representation of the used generator
    Rep& generator(Rep&) const;
	// p-adic representation of the used generator
    UTT generator() const;
	// p-adic representation of the used irreducible polynomial
    UTT irreducible() const;

	// the internal representation of the polynomial X
	// where the indeterminate is replaced by the characteristic
	// This has no meaning if exponent is 1
    Rep sage_generator() const;
    Rep indeterminate() const;
    Rep& indeterminate(Rep&) const;

	// Initialization of Elements
    Rep& init(Rep& r) const { return r=zero;}
    Rep& init( Rep&, const int32_t) const ;
    Rep& init( Rep&, const uint32_t) const ;
    Rep& init( Rep&, const int64_t) const ;
    Rep& init( Rep&, const uint64_t) const ;
    Rep& init( Rep&, const Integer) const;
    Rep& init( Rep&, const float) const ;
    Rep& init( Rep&, const double) const ;
    Rep& init( Rep& a, std::istream& s ) const { return read(a,s); }


	// Reduction of Elements
    Rep& reduce(Rep& r) const;
    Rep& reduce(Rep& r, const Rep e) const;


	// Initialization of a polynomial
    template<typename val_t, template<class,class> class Vector,template <class> class Alloc>
    Rep& init( Rep&, const Vector<val_t,Alloc<val_t> >&) const;


	// -- Misc: r <- a mod p
    Rep& assign (Rep&, const Integer) const;
    Rep& assign (Rep&, const Rep) const;
    void assign ( const size_t sz, Array r, constArray a ) const;

	// --- IO methods for the Domain
    std::istream& read ( std::istream& s );
    std::ostream& write( std::ostream& s ) const;
    std::ostream& write( std::ostream& s , const std::string& ) const;
	// --- IO methods for the Elements
    std::istream& read ( std::istream& s, Rep& a ) const;
    std::ostream& write( std::ostream& s, const Rep a ) const;

	// Conversions of the elements
    std::ostream& convert(std::ostream& s, const Rep a ) const { return write(s,a); }
    TT		convert(const Rep) const ;
    int64_t&	convert(int64_t&, const Rep) const ;
    uint64_t&	convert(uint64_t&, const Rep) const ;
    int32_t&	convert(int32_t&, const Rep) const ;
    float&	convert(float&, const Rep) const ;
    double&	convert(double&, const Rep) const ;
    uint32_t&	convert(uint32_t&, const Rep) const ;
    Integer&	convert(Integer&, const Rep) const ;

	// Test operators
	inline int operator== (const GFqDom<TT>& a) const;
	inline int operator!= (const GFqDom<TT>& a) const;

	// Miscellaneous functions
	bool areEqual(const Rep a, const Rep b) const;
	bool areNEqual(const Rep, const Rep) const;
	bool isZero(const Rep a) const;
	bool isnzero(const Rep) const;
	bool isOne(const Rep a) const;
	bool isMOne(const Rep a) const;
	bool isUnit(const Rep a) const; // Element belongs to prime subfield
	size_t length(const Rep) const;



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
	void mul (const size_t sz, Array r, constArray a, const Rep b) const;

	void div (const size_t sz, Array r, constArray a, constArray b) const;
	void div (const size_t sz, Array r, constArray a, const Rep b) const;

	void add (const size_t sz, Array r, constArray a, constArray b) const;
	void add (const size_t sz, Array r, constArray a, const Rep b) const;

	void sub (const size_t sz, Array r, constArray a, constArray b) const;
	void sub (const size_t sz, Array r, constArray a, const Rep b) const;
	void neg (const size_t sz, Array r, constArray a) const;
	void inv (const size_t sz, Array r, constArray a) const;

	Rep& axpy (Rep& r, const Rep a, const Rep b, const Rep c) const;
	void axpy (const size_t sz, Array r, const Rep a, constArray x, constArray y) const;
	void axpy (const size_t sz, Array r, const Rep a, constArray x, const Rep c) const;

	// -- axpyin: r <- r + a * x mod p
	Rep& axpyin (Rep& r, const Rep a, const Rep b) const;
	void axpyin (const size_t sz, Array r, const Rep a, constArray x) const;

	// -- axmy: r <- a * b - c mod p
	Rep& axmy (Rep& r, const Rep a, const Rep b, const Rep c) const;
	void axmy (const size_t sz, Array r, const Rep a, constArray x, constArray y) const;
	void axmy (const size_t sz, Array r, const Rep a, constArray x, const Rep c) const;

	// -- maxpy: r <- c - a * b mod p
	Rep& maxpy (Rep& r, const Rep a, const Rep b, const Rep c) const;

	// -- axmyin: r <-  a * b - r mod p
	Rep& axmyin (Rep& r, const Rep a, const Rep b) const;

	// -- maxpyin: r <- r - a * b mod p
	Rep& maxpyin (Rep& r, const Rep a, const Rep b) const;
	void maxpyin (const size_t sz, Array r, const Rep a, constArray x) const;

	// <- \sum_i a[i], return 1 if a.size() ==0,
	void reduceadd ( Rep& r, const size_t sz, constArray a ) const;

	// <- \prod_i a[i], return 1 if a.size() ==0,
	void reducemul ( Rep& r, const size_t sz, constArray a ) const;

	// <- \sum_i a[i] * b[i]
	Rep& dotprod ( Rep& r, const size_t sz, constArray a, constArray b ) const;

	// ----- random generators
	// ----- random generators

	template<class RandIter> Rep& random(RandIter& g, Rep& r) const ;
	template<class RandIter> Rep& random(RandIter& g, Rep& r, const Residu_t& s) const ;
	template<class RandIter> Rep& nonzerorandom(RandIter& g, Rep& r) const ;
	template<class RandIter> Rep& nonzerorandom(RandIter& g, Rep& r, const Residu_t& s) const ;

	typedef GIV_randIter< GFqDom<TT>, Rep> RandIter;
    typedef GeneralRingNonZeroRandIter<Self_t,RandIter> NonZeroRandIter;

#ifdef __GIVARO_COUNT__
	void clear()
	{
		_add_count = 0;
		_mul_count = 0;
		_neg_count = 0;
		_div_count = 0;
		_sub_count = 0;
		_inv_count = 0;
		_add_call = 0;
		_mul_call = 0;
		_neg_call = 0;
		_div_call = 0;
		_sub_call = 0;
		_inv_call = 0;
	}

	void info() const
	{
		std::cerr << "Mul Call: " << _mul_call << ", real: " << _mul_count << std::endl;
		std::cerr << "Add Call: " << _add_call << ", real: " << _add_count << std::endl;
		std::cerr << "Div Call: " << _div_call << ", real: " << _div_count << std::endl;
		std::cerr << "Sub Call: " << _sub_call << ", real: " << _sub_count << std::endl;
		std::cerr << "Neg Call: " << _neg_call << ", real: " << _neg_count << std::endl;
		std::cerr << "Inv Call: " << _inv_call << ", real: " << _inv_count << std::endl;
	}
#endif


#ifdef __GIVARO_COUNT__
	static    int64_t _add_count;
	static    int64_t _mul_count;
	static    int64_t _neg_count;
	static    int64_t _div_count;
	static    int64_t _sub_count;
	static    int64_t _inv_count;
	static    int64_t _add_call;
	static    int64_t _mul_call;
	static    int64_t _neg_call;
	static    int64_t _div_call;
	static    int64_t _sub_call;
	static    int64_t _inv_call;
#endif

	static void Init();
	static void End();
};

} // namespace Givaro

#include "givaro/gfq.inl"

#endif // __GIVARO_gfq1_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
