// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int_cstor.C,v $
// Copyright(c)'1994-2011 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: B. Boyer
// $Id: gmp++_rat_cstor.C,v 1.4 2011-09-16 14:28:22 bboyer Exp $
// ==========================================================================
#ifndef __GMPplusplus_rat_cstor_C__
#define __GMPplusplus_rat_cstor_C__
#include <iostream>
#include "gmp++/gmp++.h"

namespace Givaro
{
	// CONSTRUCTORS FROM INTEGERS
	Rationel::Rationel()
	{
		mpq_init((mpq_ptr)&gmp_rep);
	}

	Rationel::Rationel( const Integer & n)
	{
		mpq_init((mpq_ptr)&gmp_rep);
		mpq_set_z((mpq_ptr)&gmp_rep,(mpz_ptr)n.get_mpz());

		num = mpq_numref((mpq_ptr)&gmp_rep) ;
		den = mpq_denref((mpq_ptr)&gmp_rep) ;
	}

	Rationel::Rationel( const int  n)
	{
		mpq_init((mpq_ptr)&gmp_rep);
		mpq_set_si((mpq_ptr)&gmp_rep,(long int) n, 1L) ;

		num = mpq_numref((mpq_ptr)&gmp_rep) ;
		den = mpq_denref((mpq_ptr)&gmp_rep) ;
	}

	Rationel::Rationel( const unsigned int  n)
	{
		mpq_init((mpq_ptr)&gmp_rep);
		mpq_set_ui((mpq_ptr)&gmp_rep,(long unsigned int) n, 1UL) ;

		num = mpq_numref((mpq_ptr)&gmp_rep) ;
		den = mpq_denref((mpq_ptr)&gmp_rep) ;
	}

	Rationel::Rationel( const long int  n)
	{
		mpq_init((mpq_ptr)&gmp_rep);
		mpq_set_si((mpq_ptr)&gmp_rep, n, 1L) ;

		num = mpq_numref((mpq_ptr)&gmp_rep) ;
		den = mpq_denref((mpq_ptr)&gmp_rep) ;
	}

	Rationel::Rationel( const long unsigned int  n)
	{
		mpq_init((mpq_ptr)&gmp_rep);
		mpq_set_ui((mpq_ptr)&gmp_rep, n, 1UL) ;

		num = mpq_numref((mpq_ptr)&gmp_rep) ;
		den = mpq_denref((mpq_ptr)&gmp_rep) ;
	}

#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
	Rationel::Rationel( const long long int  n)
	{
		mpq_init((mpq_ptr)&gmp_rep);
		mpq_set_z((mpq_ptr)&gmp_rep,((Integer)n).get_mpz()) ;

		num = mpq_numref((mpq_ptr)&gmp_rep) ;
		den = mpq_denref((mpq_ptr)&gmp_rep) ;
	}

	Rationel::Rationel( const long long unsigned int  n)
	{
		mpq_init((mpq_ptr)&gmp_rep);
		mpq_set_z((mpq_ptr)&gmp_rep,((Integer)n).get_mpz()) ;

		num = mpq_numref((mpq_ptr)&gmp_rep) ;
		den = mpq_denref((mpq_ptr)&gmp_rep) ;
	}

#endif


	// CONSTRUCTORS FROM RationelS
	Rationel::Rationel( const float f, enum reduceFlag red)
	{
		mpq_init((mpq_ptr)&gmp_rep);
		mpq_set_d((mpq_ptr)&gmp_rep,(double)f);
		std::cout<< "max precision OR best approximation ?" << std::endl;
		num = mpq_numref((mpq_ptr)&gmp_rep) ;
		den = mpq_denref((mpq_ptr)&gmp_rep) ;
	}

	Rationel::Rationel( const double f, enum reduceFlag red)
	{
		mpq_init((mpq_ptr)&gmp_rep);

		mpq_set_d((mpq_ptr)&gmp_rep,f);
		// XXX
		std::cout<< "max precision OR best approximation ?" << std::endl;
		if (red)
			reduce();

		num = mpq_numref((mpq_ptr)&gmp_rep) ;
		den = mpq_denref((mpq_ptr)&gmp_rep) ;
	}

	Rationel::Rationel(const long double f, enum reduceFlag red)
	{
		mpq_init((mpq_ptr)&gmp_rep);

		// XXX
		throw "Rationel constructor from long double is not implemented yet";
		if (red)
			reduce();
		num = mpq_numref((mpq_ptr)&gmp_rep) ;
		den = mpq_denref((mpq_ptr)&gmp_rep) ;
	}

	//!@todo Initialise Rationels from mpq_f and mpfr types.

	Rationel::Rationel( const Rationel & n, enum reduceFlag red)
	{
		mpq_init((mpq_ptr)&gmp_rep);
		mpq_set( (mpq_ptr)&gmp_rep, (mpq_srcptr)n.get_mpq() );
		if (red == Reduce)
			reduce();
	}

	// CONSTRUCTORS FROM NUM AND DEN
	Rationel::Rationel( const Integer & n, const Integer & d, enum reduceFlag red)
	{
		assert(nonZero(d)); //! @todo nonZero() ?
		mpq_init((mpq_ptr)&gmp_rep);
		mpq_set_z((mpq_ptr)&gmp_rep,(mpz_ptr)n.get_mpz());
		mpq_set_den((mpq_ptr)&gmp_rep,(mpz_ptr)d.get_mpz());
		if (red == Reduce)
			reduce();

		num = mpq_numref((mpq_ptr)&gmp_rep) ;
		den = mpq_denref((mpq_ptr)&gmp_rep) ;

	}

	/*
	 *  cast | I   U   L  UL
	 *  -----------------------
	 *  I    | si  si  si ??
	 *  U    | si  ui  si ui
	 *  L    | si  si  si ??
	 *  UL   | ??  ui  ?? ui
	 */

	// castable to long int/long int
	// INT
	Rationel::Rationel( const int n, const int d,
			    enum reduceFlag red)
	{
		assert(nonZero(d));
		mpq_init((mpq_ptr)&gmp_rep);
		mpq_set_si((mpq_ptr)&gmp_rep,(long int)n, (long int)d);
		if (red == Reduce)
			reduce();

		num = mpq_numref((mpq_ptr)&gmp_rep) ;
		den = mpq_denref((mpq_ptr)&gmp_rep) ;

	}

	Rationel::Rationel( const unsigned int n, const int d,
			    enum reduceFlag red)
	{
		assert(nonZero(d));
		mpq_init((mpq_ptr)&gmp_rep);
		mpq_set_si((mpq_ptr)&gmp_rep,(long int)n, (long int)d);
		if (red == Reduce)
			reduce();

		num = mpq_numref((mpq_ptr)&gmp_rep) ;
		den = mpq_denref((mpq_ptr)&gmp_rep) ;

	}

	Rationel::Rationel( const long int n, const int d,
			    enum reduceFlag red)
	{
		assert(nonZero(d));
		mpq_init((mpq_ptr)&gmp_rep);
		mpq_set_si((mpq_ptr)&gmp_rep,n, (long int)d);
		if (red == Reduce)
			reduce();

		num = mpq_numref((mpq_ptr)&gmp_rep) ;
		den = mpq_denref((mpq_ptr)&gmp_rep) ;

	}

	Rationel::Rationel( const int n, const unsigned int d,
			    enum reduceFlag red)
	{
		assert(nonZero(d));
		mpq_init((mpq_ptr)&gmp_rep);
		mpq_set_si((mpq_ptr)&gmp_rep,(long int)n, (long int)d);
		if (red == Reduce)
			reduce();

		num = mpq_numref((mpq_ptr)&gmp_rep) ;
		den = mpq_denref((mpq_ptr)&gmp_rep) ;

	}

	Rationel::Rationel( const int n, const long int d,
			    enum reduceFlag red)
	{
		assert(nonZero(d));
		mpq_init((mpq_ptr)&gmp_rep);
		mpq_set_si((mpq_ptr)&gmp_rep,(long int)n, d);
		if (red == Reduce)
			reduce();

		num = mpq_numref((mpq_ptr)&gmp_rep) ;
		den = mpq_denref((mpq_ptr)&gmp_rep) ;

	}

	// U
	Rationel::Rationel( const unsigned int n, const long int d,
			    enum reduceFlag red)
	{
		assert(nonZero(d));
		mpq_init((mpq_ptr)&gmp_rep);
		mpq_set_si((mpq_ptr)&gmp_rep,(long int)n, d);
		if (red == Reduce)
			reduce();

		num = mpq_numref((mpq_ptr)&gmp_rep) ;
		den = mpq_denref((mpq_ptr)&gmp_rep) ;

	}
	Rationel::Rationel( const long int n, const unsigned int d,
			    enum reduceFlag red)
	{
		assert(nonZero(d));
		mpq_init((mpq_ptr)&gmp_rep);
		mpq_set_si((mpq_ptr)&gmp_rep, n, (long int)d);
		if (red == Reduce)
			reduce();

		num = mpq_numref((mpq_ptr)&gmp_rep) ;
		den = mpq_denref((mpq_ptr)&gmp_rep) ;

	}

	// L
	Rationel::Rationel( const long int n, const long int d,
			    enum reduceFlag red)
	{
		assert(nonZero(d));
		mpq_init((mpq_ptr)&gmp_rep);
		mpq_set_si((mpq_ptr)&gmp_rep, n, d);
		if (red == Reduce)
			reduce();

		num = mpq_numref((mpq_ptr)&gmp_rep) ;
		den = mpq_denref((mpq_ptr)&gmp_rep) ;

	}


	// castable to unsigned long int/unsigned long int
	// U
	Rationel::Rationel( const unsigned int n, const unsigned int d,
			    enum reduceFlag red)
	{
		assert(nonZero(d));
		mpq_init((mpq_ptr)&gmp_rep);
		mpq_set_ui((mpq_ptr)&gmp_rep,(unsigned long int)n, (unsigned long int)d);
		if (red == Reduce)
			reduce();

		num = mpq_numref((mpq_ptr)&gmp_rep) ;
		den = mpq_denref((mpq_ptr)&gmp_rep) ;

	}

	// UL
	Rationel::Rationel( const unsigned long int n, const unsigned long int d,
			    enum reduceFlag red)
	{
		assert(nonZero(d));
		mpq_init((mpq_ptr)&gmp_rep);
		mpq_set_ui((mpq_ptr)&gmp_rep,n, d);
		if (red == Reduce)
			reduce();

		num = mpq_numref((mpq_ptr)&gmp_rep) ;
		den = mpq_denref((mpq_ptr)&gmp_rep) ;

	}

	Rationel::Rationel( const unsigned long int n, const unsigned int d,
			    enum reduceFlag red)
	{
		assert(nonZero(d));
		mpq_init((mpq_ptr)&gmp_rep);
		mpq_set_ui((mpq_ptr)&gmp_rep,n, (unsigned long int)d);
		if (red == Reduce)
			reduce();

		num = mpq_numref((mpq_ptr)&gmp_rep) ;
		den = mpq_denref((mpq_ptr)&gmp_rep) ;

	}

	Rationel::Rationel( const unsigned int n, const unsigned long int d,
			    enum reduceFlag red)
	{
		assert(nonZero(d));
		mpq_init((mpq_ptr)&gmp_rep);
		mpq_set_ui((mpq_ptr)&gmp_rep,(unsigned long)n, d);
		if (red == Reduce)
			reduce();

		num = mpq_numref((mpq_ptr)&gmp_rep) ;
		den = mpq_denref((mpq_ptr)&gmp_rep) ;

	}

	// Integer fallbacks
	template<class T>
	Rationel::Rationel( const Integer & n, const T d, enum reduceFlag red)
	{
		assert(nonZero(d)); //! @todo nonZero() ?
		mpq_init((mpq_ptr)&gmp_rep);
		mpq_set_z((mpq_ptr)&gmp_rep,(mpz_ptr)n.get_mpz());
		mpq_set_den((mpq_ptr)&gmp_rep,(mpz_ptr)Integer(d).get_mpz());
		if (red == Reduce)
			reduce();

		num = mpq_numref((mpq_ptr)&gmp_rep) ;
		den = mpq_denref((mpq_ptr)&gmp_rep) ;
	}

	template<class T>
	Rationel::Rationel( const T n, const Integer & d, enum reduceFlag red)
	{
		assert(nonZero(d)); //! @todo nonZero() ?
		mpq_init((mpq_ptr)&gmp_rep);
		mpq_set_z((mpq_ptr)&gmp_rep,(mpz_ptr)Integer(n).get_mpz());
		mpq_set_den((mpq_ptr)&gmp_rep,(mpz_ptr)d.get_mpz());
		if (red == Reduce)
			reduce();

		num = mpq_numref((mpq_ptr)&gmp_rep) ;
		den = mpq_denref((mpq_ptr)&gmp_rep) ;

	}

	template<class T, class U>
	Rationel::Rationel( const T n, const U d, enum reduceFlag red)
	{
		assert(nonZero(d)); //! @todo nonZero() ?
		mpq_init((mpq_ptr)&gmp_rep);
		mpq_set_z((mpq_ptr)&gmp_rep,(mpz_ptr)Integer(n).get_mpz());
		mpq_set_den((mpq_ptr)&gmp_rep,(mpz_ptr)Integer(d).get_mpz());
		if (red == Reduce)
			reduce();

		num = mpq_numref((mpq_ptr)&gmp_rep) ;
		den = mpq_denref((mpq_ptr)&gmp_rep) ;

	}


}
#endif // __GMPplusplus_rat_cstor_C__

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen
