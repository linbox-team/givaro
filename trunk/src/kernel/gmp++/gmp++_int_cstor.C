// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int_cstor.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: M. Samama, T. Gautier
// $Id: gmp++_int_cstor.C,v 1.4 2009-09-17 14:28:22 jgdumas Exp $
// ==========================================================================
#ifndef __GMPplusplus_CSTOR_C__
#define __GMPplusplus_CSTOR_C__
#include <iostream>
#include "gmp++/gmp++.h"


namespace Givaro {
	//------------------------------------- predefined null and one
	const Integer Integer::zero(0UL);
	const Integer Integer::one(1UL);


	// -- Integer(const char *s)
	Integer::Integer(const char *s)
	{
		mpz_init_set_str((mpz_ptr)&gmp_rep, s, 10);
	}


	Integer& Integer::copy(const Integer &n)
	{
		if (this == &n) return *this;
		mpz_set ( (mpz_ptr)&gmp_rep, (mpz_srcptr)&(n.gmp_rep)) ;
		return *this ;
	}

	void importWords(Integer& x, size_t count, int order, int size, int endian, size_t nails, const void* op)
	{
		mpz_import( (mpz_ptr)&(x.gmp_rep), count, order, size, endian, nails, op);
	}

}

#endif

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen
