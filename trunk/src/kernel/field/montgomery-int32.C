// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givmontg32.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givmontg32.C,v 1.6 2011-02-02 17:16:43 bboyer Exp $
// ==========================================================================
// Description:

#ifndef __GIVARO_zpz_givmontg32_C
#define __GIVARO_zpz_givmontg32_C

#include <iostream>
#include "givaro/montgomery-int32.h"

namespace Givaro {
  // Returns d, and u and v such that u a + v b = d;
// JGD 04.11.1999
// const int32_t Montgomery<int32_t>::gcdext
//   ( int32_t& u, int32_t& v, const int32_t a, const int32_t b )
int32_t& Montgomery<int32_t>::gcdext
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
   return d=int32_t(u3);
//    return u3;
}

int32_t& Montgomery<int32_t>::invext
  ( int32_t& u, const int32_t a, const int32_t b ) const
{
    long u1,u3;
    long v1,v3;
   u1 = 1;
   u3 = long(a);
   v1 = 0;
   v3 = long(b);
   while (v3 != 0)
   {
	 long q, t1, t3;
	q = u3 / v3;
	t1 = u1 - q * v1; t3 = u3 - q * v3;
	u1 = v1; u3 = v3; v1 = t1; v3 = t3;
   }
   if (u1 < 0)
       return u = int32_t(u1+b);
   else
       return u = int32_t(u1);
}

int32_t Montgomery<int32_t>::invext
  ( const int32_t a, const int32_t b ) const
{
    int32_t tmp;
    return invext(tmp, a, b);
}


void Montgomery<int32_t>::Init()
{
}

void Montgomery<int32_t>::End()
{
}
} // namespace Givaro

#endif
