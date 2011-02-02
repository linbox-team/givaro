// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givmontg32.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givmontg32.C,v 1.5 2011-02-02 16:23:56 bboyer Exp $
// ==========================================================================
// Description:

#include <iostream>
#include "givaro/givmontg32.h"

  // Returns d, and u and v such that u a + v b = d;
// JGD 04.11.1999
// const int32 Montgomery<Std32>::gcdext
//   ( int32& u, int32& v, const int32 a, const int32 b )
int32& Montgomery<Std32>::gcdext
  ( int32& d,  int32& u, int32& v, const int32 a, const int32 b ) const
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
   u = u1;
   v = u2;
   return d=u3;
//    return u3;
}

int32& Montgomery<Std32>::invext
  ( int32& u, const int32 a, const int32 b ) const
{
    long u1,u3;
    long v1,v3;
   u1 = 1; u3 = a;
   v1 = 0; v3 = b;
   while (v3 != 0)
   {
	 long q, t1, t3;
	q = u3 / v3;
	t1 = u1 - q * v1; t3 = u3 - q * v3;
	u1 = v1; u3 = v3; v1 = t1; v3 = t3;
   }
   if (u1 < 0)
       return u = u1+b;
   else
       return u = u1;
}

int32 Montgomery<Std32>::invext
  ( const int32 a, const int32 b ) const
{
    int32 tmp;
    return invext(tmp, a, b);
}


void Montgomery<Std32>::Init()
{
}

void Montgomery<Std32>::End()
{
}

