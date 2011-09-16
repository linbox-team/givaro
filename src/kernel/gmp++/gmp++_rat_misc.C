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
#ifndef __GMPplusplus_rat_misc_C__
#define __GMPplusplus_rat_misc_C__
#include <iostream>
#include "gmp++/gmp++.h"

namespace Givaro
{
	Integer  Rationel::getDenom() const
	{
		Integer d ;
		mpq_get_den(d.get_mpz(), (mpq_srcptr)&gmp_rep);
		return d;
	}

	Integer  Rationel::getNumer() const
	{
		Integer n ;
		mpq_get_num(n.get_mpz(), (mpq_srcptr)&gmp_rep);
		return n;
	}

	mpq_ptr Rationel::get_mpq() const
	{
		return (mpq_ptr)&gmp_rep;
	}

	mpz_ptr Rationel::get_mpq_den() const
	{
		return (mpz_ptr)den;
	}

	mpz_ptr Rationel::get_mpq_num() const
	{
		return (mpz_ptr)num;
	}

	Rationel& Rationel::reduce()
	{
		mpq_canonicalize( (mpq_ptr)&gmp_rep );
		return *this ;
	}

	Rationel& Rationel::reduce(Rationel & r) //const
	{
		mpq_canonicalize( (mpq_ptr)(r.get_mpq()) );
		return r ;
	}

}

#endif // __GMPplusplus_rat_misc_C__

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen
