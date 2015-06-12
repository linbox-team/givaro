// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/rational/givratreconstruct.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: Jean-Guillaume Dumas
// $Id: givratreconstruct.C,v 1.4 2010-04-12 15:54:39 jgdumas Exp $
// ==========================================================================
// Description:

#include "givaro/givrational.h"
#include <iostream>

namespace Givaro {


	Rational::Rational(const Integer& f, const Integer& m, const Integer& k, bool recurs )
	{
		bool res = this->ratrecon(f,m,k,recurs);
		if (recurs)
			for( Integer newk = k + 1L; (!res) && (newk<f) ; ++newk)
				res = this->ratrecon(f,m,newk,true);
	}



	bool Rational::ratrecon(const Integer& f, const Integer& m, const Integer& k, bool recurs )
	{

#ifdef DEBUG
		std::cerr << "RatRecon : " << f << " " << m << " " << k << std::endl;
#endif


		Integer r0, t0, r1, t1, q, u;
		r0=m;
		t0=0;
		r1=f;
		if (f<0) r1+= m;
		t1=1;
		while(r1>=k)
		{
#ifdef DEBUG
			std::cerr << "r0: " << r0
			<< ", r1: " << r1
			<< ", q: " << q
			<< ", t0: " << t0
			<< ", t1: " << t1
			<< std::endl;
#endif
			q = r0;
			q /= r1;        // r0/r1

			u = r1;
			r1 = r0;  	// r1 <-- r0
			r0 = u;	        // r0 <-- r1
			u *= q;
			r1 -= u;	// r1 <-- r0-q*r1
			if (r1 == 0) break;

			u = t1;
			t1 = t0;  	// r1 <-- r0
			t0 = u;	        // r0 <-- r1
			u *= q;
			t1 -= u;	// r1 <-- r0-q*r1
		}

		// [GG, MCA, 1999] Theorem 5.26
		// (i)
		if (t1 < 0) {
			num = -r1;
			den = -t1;
		} else {
			num = r1;
			den = t1;
		}

		// (ii)
		if (gcd(num,den) != 1) {
			//         std::cerr << "r0: " << r0
			//                   << ", r1: " << r1
			//                   << ", k: " << k
			//                   << std::endl;

			Integer gar1, gar2;
			for( q = 1, gar1 = r0-r1, gar2 = r0 ; (gar1 >= k) || (gar2<k); ++q ) {
				gar1 -= r1;
				gar2 -= r1;
			}

			r0 -= q * r1;
			t0 -= q * t1;

			if (t0 < 0) {
				num = -r0;
				den = -t0;
			} else {
				num = r0;
				den = t0;
			}

			if (t0 > m/k) {
				if (!recurs)
					std::cerr
					<< "*** Error *** No rational reconstruction of "
					<< f
					<< " modulo "
					<< m
					<< " with denominator <= "
					<< (m/k)
					<< std::endl;
			}
			if (gcd(num,den) != 1) {
				if (!recurs)
					std::cerr
					<< "*** Error *** There exists no rational reconstruction of "
					<< f
					<< " modulo "
					<< m
					<< " with |numerator| < "
					<< k
					<< std::endl
					<< "*** Error *** But "
					<< num
					<< " = "
					<< den
					<< " * "
					<< f
					<< " modulo "
					<< m
					<< std::endl;
				return false;
			}
		}
		// std::cerr << "RatRecon End " << std::endl;
		return true;
	}


} // namespace Givaro

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen
