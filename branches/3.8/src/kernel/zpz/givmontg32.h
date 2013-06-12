// ==========================================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// author: JG Dumas (from P. Zimmermann's Montgomery implementation)
// $Id: givmontg32.h,v 1.15 2011-02-04 14:11:46 jgdumas Exp $
// ==========================================================================

/*! @file givmontg32.h
 * @ingroup zpz
 * @brief NO DOC
 */

#ifndef __GIVARO_montg32_H
#define __GIVARO_montg32_H

#include "givaro/givbasictype.h"
#include "givaro/giverror.h"
#include "givaro/giv_randiter.h"
#include <math.h>



#define B32 65536UL
#define MASK32 65535UL
#define HALF_BITS32 16

namespace Givaro {

	template<class TYPE> class Montgomery;


	/*! @brief This class implements the standard arithmetic with Modulo Elements.
	 *   Reduction is made through Montgomery's reduction.
	 *   Representation of a is by storing (aB).
	 *   - We must have gcd(p,2)=1
	 *   - We must have \f$(p-1)^2 + p(B-1) < B^2 \f$, i.e. \f$2<p \leq 40504\f$ for \f$B=2^16\f$.
     *   - m max is 40503
     *   - p max is 40499
	 */
	template<>
	class Montgomery<Std32> {
	public:
		// ----- Exported Types and constantes
		typedef uint32_t Residu_t;                    // - type to store residue
		enum { size_rep = sizeof(Residu_t) };      // - size of the storage type
		// ----- Representation of Element of the domain Montgomery
		typedef uint32_t Rep;
		typedef uint32_t Element;

		// ----- Constructor
		Montgomery() : _p(0UL), _dp(0.0), zero(0UL), one(1UL), mOne(0UL) {}

		Montgomery( Residu_t p, int = 1) :
			_p(  (Residu_t)  p),
			_Bp( (Residu_t)  B32%p),
			_B2p((Residu_t)  (_Bp<<HALF_BITS32) % p),
			_B3p((Residu_t)  (_B2p<<HALF_BITS32) % p),
			_nim((Residu_t)  -Montgomery<Std32>::invext((int32_t)_p,B32) ),
			_dp( (double)    p),
			zero((Residu_t)  0UL),
			one( (Residu_t)  redcsal( _B2p ) ),
			mOne( _p - one )
		{}

		Montgomery( const Montgomery<Std32>& F)
		: _p(F._p), _Bp(F._Bp), _B2p( F._B2p), _B3p( F._B3p), _nim(F._nim),_dp(F._dp), zero(0UL), one(F.one),mOne(F.mOne)
		{ }


		int operator==( const Montgomery<Std32>& BC) const
		{ return _p == BC._p;}
		int operator!=( const Montgomery<Std32>& BC) const
		{ return _p != BC._p;}

		Montgomery<Std32>& operator=( const Montgomery<Std32>& F)
		{
			this->_p = F._p;
			this->_Bp = F._Bp;
			this->_B2p = F._B2p;
			this->_B3p = F._B3p;
			this->_nim = F._nim;
			this->_dp = F._dp;
			assign(const_cast<Element&>(one),F.one);
			assign(const_cast<Element&>(mOne),F.mOne);
			assign(const_cast<Element&>(zero),F.zero);
			return *this;
		}

		// ----- Access to the modulus
		Residu_t residu() const;
		Residu_t size() const
		{return _p;}
		Rep access( const Rep a ) const
		{ return a; }
		Residu_t characteristic() const
		{ return _p; }
		Residu_t characteristic(Residu_t p) const
		{ return p=_p; }
		Residu_t cardinality() const
		{ return _p; }


		// ----- Access to the modulus
		Rep& init( Rep& a ) const;
		Rep& init( Rep& r, const long a) const ;
		Rep& init( Rep& r, const unsigned long a) const ;
		Rep& init( Rep& a, const int i) const ;
		Rep& init( Rep& a, const unsigned int i) const ;
		Rep& init ( Rep& r, const Integer& residu ) const ;

		// Initialisation from double ( added for FFLAS usage) (C Pernet)
		Rep& init( Rep& a, const double i) const;
		Rep& init( Rep& a, const float i) const;

		unsigned long int& convert(unsigned long int& r, const Rep a) const
		{
			uint32_t ur;
			return r = (unsigned long)redc(ur,a);}

			uint32_t& convert(uint32_t& r, const Rep a) const
			{
				unsigned long ur;
				return r = (uint32_t)convert(ur, a);
			}

			int32_t& convert(int32_t& r, const Rep a) const
			{
				unsigned long ur;
				return r = (int32_t)convert(ur, a);
			}

			long int& convert(long int& r, const Rep a) const
			{
				unsigned long ur;
				return r = (long int)convert(ur, a);
			}

			Integer& convert(Integer& i, const Rep a) const
			{
				unsigned long ur;
				return i = (Integer)convert(ur, a);
			}

			// Conversion to double ( added for FFLAS usage) (C Pernet)
			float& convert(float& r, const Rep a ) const
			{
				unsigned long ur;
				return r = (float)convert(ur, a); }
				double& convert(double& r, const Rep a ) const
				{
					unsigned long ur;
					return r = (double)convert(ur, a); }

					// ----- Misc methods
					int isZero( const Rep a ) const;
					int isOne ( const Rep a ) const;
					size_t length ( const Rep a ) const;

					// ----- Equality between two Elements
					int areEqual(const  Rep& a, const Rep& b) const
					{
						return a==b;
					}

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

					// -- axpy: r <- a * x + y mod p
					Rep& axpy  (Rep& r, const Rep a, const Rep b, const Rep c) const;
					// -- axpyin: r <- r + a * x mod p
					Rep& axpyin(Rep& r, const Rep a, const Rep b) const;
					// -- axmy: r <- a * x - y mod p
					Rep& axmy  (Rep& r, const Rep a, const Rep b, const Rep c) const;
					// -- axmyin: r <- a * x - r  mod p
					Rep& axmyin(Rep& r, const Rep a, const Rep b) const;
					// -- maxpy: r <- c - a * b mod p
					Rep& maxpy  (Rep& r, const Rep a, const Rep b, const Rep c) const;
					// -- maxpyin: r <- r - a * x mod p
					Rep& maxpyin(Rep& r, const Rep a, const Rep b) const;
					// -- Misc: r <- a mod p
					Rep& assign ( Rep& r, const Rep a) const;

					// ----- random generators
					template< class RandIter > Rep& random(RandIter&, Rep& r) const ;
					template< class RandIter > Rep& random(RandIter&, Rep& r, long s) const ;
					template< class RandIter > Rep& random(RandIter&, Rep& r, const Rep& b) const ;
					template< class RandIter > Rep& nonzerorandom(RandIter&, Rep& r) const ;
					template< class RandIter > Rep& nonzerorandom(RandIter&, Rep& r, long s) const ;
					template< class RandIter > Rep& nonzerorandom(RandIter&, Rep& r, const Rep& b) const ;

					typedef GIV_randIter< Montgomery<Std32> , Rep > randIter;

					// --- IO methods
					std::istream& read ( std::istream& s );
					std::ostream& write( std::ostream& s ) const;
					std::istream& read ( std::istream& s, Rep& a ) const;
					std::ostream& write( std::ostream& s, const Rep a ) const;

	protected:
					// -- based on modular inverse, d = a*u + b*v
					//   static const int32_t gcdext ( int32_t& u, int32_t& v, const int32_t a, const int32_t b );
					int32_t& gcdext (int32_t& d, int32_t& u, int32_t& v, const int32_t a, const int32_t b ) const;
					int32_t& invext (int32_t& u, const int32_t a, const int32_t b ) const;
					int32_t invext(const int32_t a, const int32_t b ) const;



					Element& redc(Element&, const Element) const ;
					Element redcal(const Element) const;
					Element redcsal(const Element) const;
					Element& redcin(Element&) const;
					Element& redcs(Element&, const Element) const;
					Element& redcsin(Element&) const;


	protected:
					// -- data representation of the domain:
					Residu_t _p;
					Residu_t _Bp;
					Residu_t _B2p;
					Residu_t _B3p;
					Residu_t _nim;
					double _dp;


					static void Init();
					static void End();

	public:
					// ----- Constantes
					const Rep zero;
					const Rep one;
					const Rep mOne;

    public: static inline Residu_t getMaxModulus() { return 40503; }
    
	};


} // namespace Givaro

#include "givaro/givmontg32.inl"

#endif //  __GIVARO_montg32_H

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
