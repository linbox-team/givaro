// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givzpz32uns.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givzpz32uns.C,v 1.5 2011-02-02 17:16:43 bboyer Exp $
// ==========================================================================
// Description:

#include <iostream>
#include "givaro/givzpz32uns.h"

namespace Givaro {

	// Returns d, and u and v such that u a + v b = d;
	// JGD 04.11.1999
	// const int32_t ZpzDom<Unsigned32>::gcdext
	//   ( int32_t& u, int32_t& v, const int32_t a, const int32_t b )
	int32_t& ZpzDom<Unsigned32>::gcdext
	( int32_t& d,  int32_t& u, int32_t& v, const int32_t a, const int32_t b ) const
	{
		long u1,u2,u3;
		long v1,v2,v3;
		u1 = 1; u2 = 0; u3 = a;
		v1 = 0; v2 = 1; v3 = b;
		while (v3 != 0)
		{
			long q , t1, t2 ,t3;
			q = u3 / v3;
			t1 = u1 - q * v1; t2 = u2 - q * v2; t3 = u3 - q * v3;
			u1 = v1; u2 = v2; u3 = v3; v1 = t1; v2 = t2; v3 = t3;
		}
		u = int32_t(u1);
		v = int32_t(u2);
		return d= int32_t(u3);
		//    return u3;
	}

	uint32_t& ZpzDom<Unsigned32>::invext
	( uint32_t& u1, const uint32_t a, const uint32_t b ) const
	{
            u1=one;
            uint32_t r0(_p), r1(a);
            uint32_t q(r0/r1);
            
            r0 -= q * r1;
            if (r0 == zero) return u1;
            uint32_t u0 = q;
            
            q = r1/r0;
            r1 -= q * r0; 
            
            while (r1 != zero) {
                u1 += q * u0;
                
                q = r0/r1;
                r0 -= q * r1;
                if (r0 == zero) return u1;
                u0 += q * u1;
                
                q = r1/r0;
                r1 -= q * r0; 
                
            }
            
            return u1=_p-u0;
	}


	void ZpzDom<Unsigned32>::Init()
	{
	}

	void ZpzDom<Unsigned32>::End()
	{
	}

} // namespace Givaro

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
