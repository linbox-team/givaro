// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int_io.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors:  B. Boyer
// $Id: gmp++_int_io.C,v 1.7 2011-09-17 14:28:22 bboyer Exp $
// ==========================================================================
// Description:

#ifndef __GIVARO_gmpxx_gmpxx_rat_io_C
#define __GIVARO_gmpxx_gmpxx_rat_io_C
#include <iostream>
#include <stdlib.h>
#include <sstream>
#ifndef __GIVARO_INLINE_ALL
#include "gmp++/gmp++.h"
#endif

namespace Givaro {

	// Sortie nonsignee : 321321 meme si n = -321321, par exemple
	std::ostream& absOutput(std::ostream &o, const Rationel &n)
	{
		mpq_out_str(reinterpret_cast<FILE*>(&o),10,(mpq_srcptr)abs(n).get_mpq());
		return o;
	}

	// Sortie signee : +321321 ou -321321, par exemple
	std::ostream& Rationel::print(std::ostream &o) const
	{
		// mpq_out_str(reinterpret_cast<FILE*>(&o),10,(mpq_srcptr)&gmp_rep);
		// return o;

		return o << (mpq_srcptr)&gmp_rep;
	}

	Rationel::operator std::string () const
	{
		std::ostringstream s;
		print(s);
		return s.str();
	}

#if 0
	Rationel::Rationel(const std::vector<mp_limb_t>& v)
	{
		size_t s = v.size();
		if (s) {
			mpz_init_set_ui((mpz_ptr)&gmp_rep, v[0]);
			Integer base(256), prod, tmp;
			prod = base = pow(base, (unsigned long)sizeof(mp_limb_t) );

			std::vector<mp_limb_t>::const_iterator vi = v.begin();
			for(++vi;vi != v.end();++vi) {
				mpz_mul_ui( (mpz_ptr)&tmp.gmp_rep, (mpz_ptr)&prod.gmp_rep, *vi);
				*this += tmp;
				prod *= base;
			}
		} else
			mpz_init( (mpz_ptr)&gmp_rep );

	}

	Rationel::operator std::vector<mp_limb_t> () const
	{
		size_t s = mpz_size( (mpz_srcptr)&(gmp_rep) );
		std::vector<mp_limb_t> v(s);
		std::vector<mp_limb_t>::iterator vi = v.begin();
		for(mp_size_t i = 0;vi != v.end();++vi, ++i)
			*vi = mpz_getlimbn( (mpz_srcptr)& (gmp_rep) ,i);
		return v;
	}
#endif

	// Entree au format de la sortie
	std::istream& operator>> (std::istream& inp, Rationel& a)
	{
		return inp >>  (mpq_ptr)a.get_mpq();
	}

	std::ostream& operator<< (std::ostream& outp, const Rationel& a)
	{
		return a.print(outp) ;
	}

}

#endif // __GIVARO_gmpxx_gmpxx_rat_io_C

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen
