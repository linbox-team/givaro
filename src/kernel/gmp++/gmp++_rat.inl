// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int_cstor.C,v $
// Copyright(c)'1994-2011 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: B. Boyer
// $Id: gmp++_rat_cstor.C,v 1.4 2011-09-16 14:28:22 briceboyer Exp $
// ==========================================================================

#ifndef __GIVARO_gmpxx_gmpxx_rat_INL
#define __GIVARO_gmpxx_gmpxx_rat_INL

namespace Givaro
{
	// Integer constructors fallbacks
	template<class T>
	Rationel::Rationel( Integer & n, T d, enum reduceFlag red)
	{
		assert(nonZero(d)); //! @todo nonZero() ?
		mpq_init((mpq_ptr)&gmp_rep);
		mpq_set_z((mpq_ptr)&gmp_rep,(mpz_srcptr)n.get_mpz_const());
		mpq_set_den((mpq_ptr)&gmp_rep,(mpz_ptr)Integer(d).get_mpz());
		if (red == Reduce)
			reduce();

		num = mpq_numref((mpq_ptr)&gmp_rep) ;
		den = mpq_denref((mpq_ptr)&gmp_rep) ;
	}

	template<class T>
	Rationel::Rationel( T n, Integer & d, enum reduceFlag red)
	{
		assert(nonZero(d)); //! @todo nonZero() ?
		mpq_init((mpq_ptr)&gmp_rep);
		mpq_set_z((mpq_ptr)&gmp_rep,(mpz_ptr)Integer(n).get_mpz());
		mpq_set_den((mpq_ptr)&gmp_rep,(mpz_srcptr)d.get_mpz_const());
		if (red == Reduce)
			reduce();

		num = mpq_numref((mpq_ptr)&gmp_rep) ;
		den = mpq_denref((mpq_ptr)&gmp_rep) ;

	}

	template<class T, class U>
	Rationel::Rationel( T n, U d, enum reduceFlag red)
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

	template<class T, class U>
	Rationel::Rationel( T & n, U & d, enum reduceFlag red)
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


#endif // __GIVARO_gmpxx_gmpxx_rat_INL
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen
