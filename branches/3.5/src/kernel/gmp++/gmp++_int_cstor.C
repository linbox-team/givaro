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
#ifndef __GIVARO_INLINE_ALL
#include "gmp++/gmp++.h"
#endif


namespace Givaro {
	//-----------------------------~Integer()
	Integer::~Integer()
	{
		mpz_clear((mpz_ptr)&gmp_rep) ;
	}

	//-------------------------------Integer(const Integer &n)
	Integer::Integer(const Integer &n)
	{
		mpz_init_set ( (mpz_ptr)&gmp_rep, (mpz_srcptr)&(n.gmp_rep)) ;
	}

	//-----------------------------Integer(int n)
	Integer::Integer(int n)
	{
		mpz_init_set_si((mpz_ptr)&gmp_rep, n) ;
	}

	//-----------------------------Integer(uint n)
	Integer::Integer(unsigned char n)
	{
		mpz_init_set_ui((mpz_ptr)&gmp_rep, n) ;
	}

	//-----------------------------Integer(uint n)
	Integer::Integer(unsigned int n)
	{
		mpz_init_set_ui((mpz_ptr)&gmp_rep, n) ;
	}

	//-----------------------------Integer(long n)
	Integer::Integer(long n)
	{
		mpz_init_set_si((mpz_ptr)&gmp_rep, n) ;
	}

	//-----------------------------Integer(unsigned long n)
	Integer::Integer(unsigned long n)
	{
		mpz_init_set_ui((mpz_ptr)&gmp_rep, n) ;
	}

#ifdef __USE_GMPPLUSPLUS_SIXTYFOUR__
#include <stdio.h>
	//-----------------------------Integer(long long n)
	// log[10](2^8) < 2.408239966
	Integer::Integer(long long n)
	{
		char * tmp = new char[long(2.408239966*sizeof(long long))+1]; sprintf(tmp,"%lld",n);
		mpz_init_set_str((mpz_ptr)&gmp_rep, tmp, 10) ;
		delete [] tmp;
	}

	//-----------------------------Integer(unsigned long long n)
	// log[10](2^8) < 2.408239966
	Integer::Integer(unsigned long long n)
	{
		char * tmp = new char[ long(2.408239966*sizeof(unsigned long long))+1];
		sprintf(tmp,"%llu",n);
		mpz_init_set_str((mpz_ptr)&gmp_rep, tmp, 10) ;
		delete [] tmp;
	}
#endif


	//-----------------------------Integer(double)
	Integer::Integer(double d)
	{
		mpz_init_set_d((mpz_ptr)&gmp_rep, d) ;
	}

	// -- Integer(const char *s)
	Integer::Integer(const char *s)
	{
		mpz_init_set_str((mpz_ptr)&gmp_rep, s, 10);
	}

	//------------------------------------------operator = (const Integer &n)
	Integer& Integer::logcpy(const Integer &n)
	{
		if (this == &n) return *this;
		mpz_set ( (mpz_ptr)&gmp_rep, (mpz_srcptr)&(n.gmp_rep)) ;
		return *this;
	}

	// same as logcopy
	Integer& Integer::operator = (const Integer &n)
	{
		return logcpy(n) ;
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
