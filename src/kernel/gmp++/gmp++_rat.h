// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_Rationel.h,v $
// Copyright(c)'1994-2011 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: B. Boyer <bboyer@imag.fr>
// $Id: givRationel.h,v 0.0 2011-09-15 18:23:56 bboyer Exp $
// ==========================================================================

#ifndef __GIVARO_GMPplusplus_Rationel_H
#define __GIVARO_GMPplusplus_Rationel_H
// #define __GIVARO_Rationel_H

/*! @file kernel/gmp++/gmp++_Rationel.h
 * @ingroup Rationel
 * Core Rationel from GMP
 */

#include <gmp++/gmp++.h>

namespace Givaro
{
	class Rationel ;
	Rationel 	abs(const Rationel& n);
	std::istream& 	operator >> (std::istream &i, Rationel& n);
	std::ostream& 	operator << (std::ostream &o, const Rationel& n);
	std::ostream& 	absOutput (std::ostream &o, const Rationel& n);
	void 		importWords(Rationel&, size_t, int, int, int, size_t, const void*);


	class Rationel
	{
		public:
			enum reduceFlag  { Reduce = 0x1, NoReduce = 0x0 } ;
			// enum ReducedFlag { Reduced = 0x2, NoReduced = 0x3 } ;

		protected:

			typedef __mpq_struct Rep;
			typedef __mpz_struct RawRep ;

			Rep gmp_rep;
			RawRep * num ;
			RawRep * den ;

			int privSign() const;

			const Rep* get_rep() const
			{
				return &gmp_rep;
			}

		public :

			// CONSTRUCTORS (gmp++_rat_cstor.C)
			//! default constructor (0/1).
			Rationel() ;

			//! constructors from a single numerator.
			//!@param n numerator, will make fraction n/1
			//@{

			Rationel( Integer & n) ;

			Rationel( int  n) ;
			Rationel( unsigned int  n) ;

			Rationel( long int  n) ;
			Rationel( long unsigned int  n) ;

#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
			Rationel( long long int  n) ;
			Rationel( long long unsigned int  n) ;
#endif
			//@}

					//! constructors from a  numerator and a denominator.
			/*! @param n numerator
			 * @param d denominator, will produce fraction \c n/d
			 * @param reduceFlag by defaut, noting is reduced.
			 * @pre we suppose \c d!=0. this is not asserted in \c NDEBUG mode !
			 */
			//@{
			Rationel( Integer & n, Integer & d,
				  enum reduceFlag = NoReduce) ;
			Rationel( int n, int d,
				  enum reduceFlag = NoReduce) ;
			Rationel( unsigned int n, int d,
				  enum reduceFlag = NoReduce) ;
			Rationel( long int n, int d,
				  enum reduceFlag = NoReduce) ;
			Rationel( long unsigned int n, int d,
				  enum reduceFlag = NoReduce) ;

			Rationel( int n, unsigned int d,
				  enum reduceFlag = NoReduce) ;
			Rationel( int n, long int d,
				  enum reduceFlag = NoReduce) ;
			Rationel( int n, unsigned long int d,
				  enum reduceFlag = NoReduce) ;

			Rationel( unsigned int n, unsigned int d,
				  enum reduceFlag = NoReduce);
			Rationel( unsigned int n, long int d,
				  enum reduceFlag = NoReduce) ;
			Rationel( unsigned int n, unsigned long int d,
				  enum reduceFlag = NoReduce);

			Rationel( unsigned long int n, unsigned int d,
				  enum reduceFlag = NoReduce) ;
			Rationel( long int n, unsigned int d,
				  enum reduceFlag = NoReduce) ;


			Rationel( long int n, long int d,
				  enum reduceFlag = NoReduce) ;
			Rationel( unsigned long int n, long int d,
				  enum reduceFlag = NoReduce) ;
			Rationel( long int n, unsigned long int d,
				  enum reduceFlag = NoReduce) ;

			Rationel( unsigned long int n, unsigned long int d,
				  enum reduceFlag = NoReduce) ;

			template<class T>
			Rationel( Integer & n, T d,
				  enum reduceFlag = NoReduce) ;
			template<class T>
			Rationel( T n, Integer & d,
				  enum reduceFlag = NoReduce) ;
			template<class T, class U>
			Rationel( T n, U d, enum reduceFlag = NoReduce) ;
			template<class T, class U>
			Rationel( T & n, U & d, enum reduceFlag = NoReduce) ;

			//@}

			//! constructors from another Rationel
			/*! @param f the Rationel to be represented
			 * @param reduceFlag a flag to start reduction or not.
			 */
			//@{
			Rationel( Rationel & f, enum reduceFlag = NoReduce) ;

			Rationel( float f, enum reduceFlag = NoReduce) ;
			Rationel( double f, enum reduceFlag = NoReduce) ;
			Rationel( long double f, enum reduceFlag = NoReduce) ;
			//@}

			// destructors
			/*! Destructor.
			 * Clearing the \c mpq representation.
			 */
			~Rationel()
			{
				mpq_clear((mpq_ptr)&gmp_rep);
			}

			// ADD

			static Rationel& addin(Rationel& res, const Rationel& n);
			static Rationel& addin(Rationel& res, const Integer& n) ;

			// SUB
			static Rationel& negin (Rationel& res)
			{
				mpq_neg((mpq_ptr)res.get_mpq(),
					(mpq_ptr)res.get_mpq());
				return res ;
			}
			Rationel& negin ()
			{
				mpq_neg((mpq_ptr)&gmp_rep,(mpq_ptr)&gmp_rep);
				return *this;
			}


			// Conversions/Casts (gmp++_rat_cast.C)
			operator std::string() const ;

			// (gmp++_rat_compare.C)
			friend Rationel abs(const Rationel& n);

			// Input/output of Rationels (gmp++_rat_io.C)
			/*! @name I/O
			*/
			//@{
			friend std::istream& operator >> (std::istream &i, Rationel & n);
			friend std::ostream& operator << (std::ostream &o, const Rationel & n);
			friend std::ostream& absOutput (std::ostream &o, const Rationel& n);

			// friend void importWords(Integer&, size_t, int, int, int, size_t, const void*);

			std::ostream& print( std::ostream& o ) const;
			//@}


			// basic acessing operations and other tools (gmp++_rat_misc.C)
			/*! Gets the denominator of a \c Rationel.
			 * @return the \c Integer denominator.
			 */
			Integer  getDenom() const ;
			/*! Gets the numerator of a \c Rationel.
			 * @return the \c Integer numerator .
			 */
			Integer  getNumer() const ;

			/*! Retrieve the GMP representation of a Rationel.
			 * @return a pointer to this representation.
			 */
			mpq_ptr get_mpq()     const ;
			/*! Retrieve the GMP representation of the denominator of a Rationel.
			 * @return a pointer to this denominator (integer).
			 */
			mpz_ptr get_mpq_den() const ;
			/*! Retrieve the GMP representation of the numerator of a Rationel.
			 * @return a pointer to this numerator (integer).
			 */
			mpz_ptr get_mpq_num() const ;

			/*! Reduces (inplace) a fraction to a canonical representation.
			 * @return a reference to self.
			 */
			Rationel& reduce();
			/*! Reduces a fraction to a canonical representation.
			 * @param r a \c Rationel
			 * @return a reference to the reduced \p r.
			 */
			static Rationel& reduce(Rationel & r) ;//const

			static inline int isZero(const Rationel &n)
			{
				return (mpq_sgn((mpq_srcptr)&n.gmp_rep) == 0) ;
			}

			int isZero()
			{
				return (mpq_sgn((mpq_srcptr)&gmp_rep) == 0);
			}

			//------------------------------------------operator = (const Integer &n)
			Rationel& logcpy(const Rationel &n)
			{
				if (this == &n) return *this;
				mpq_set ( (mpq_ptr)&gmp_rep, (mpq_srcptr)&(n.gmp_rep)) ;
				return *this;
			}

			// same as logcopy
			Rationel& operator = (const Rationel &n)
			{
				return logcpy(n) ;
			}

			Rationel& operator = (const Integer &n)
			{
				mpq_set_z ( (mpq_ptr)&gmp_rep, (mpz_srcptr)(n.get_mpz_const())) ;
				return *this;
			}

			static void setInteger(Rationel &f, const Integer & n)
			{
				// (( (mpq_ptr)&f.gmp_rep )->_mp_num ) =  (mpz_srcptr) n.get_mpz_const() ;
				// mpz_t a ;
				// mpz_init(a);
				mpz_ptr a = mpq_numref( (mpq_ptr)&f.gmp_rep );
				a = (mpz_ptr) n.get_mpz_const() ;
			}

		protected:
			static reduceFlag flags ;    //!< flag that indicates reduction is done or not after an operation.   By default, this is Reduce (as in GMP).



	};
}

#include "gmp++_rat.inl"

#endif // __GIVARO_GMPplusplus_Rationel_H

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen
