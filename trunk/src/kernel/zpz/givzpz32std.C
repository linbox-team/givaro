// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givzpz32std.C,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givzpz32std.C,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description:

#include <iostream>
#include "givaro/givzpz32std.h"

  // Returns d, and u and v such that u a + v b = d;
// JGD 04.11.1999
// const int32 ZpzDom<Std32>::gcdext 
//   ( int32& u, int32& v, const int32 a, const int32 b )
int32& ZpzDom<Std32>::gcdext 
  ( int32& d,  int32& u, int32& v, const int32 a, const int32 b ) const
{
   register long u1,u2,u3;
   register long v1,v2,v3;
   u1 = 1; u2 = 0; u3 = a;
   v1 = 0; v2 = 1; v3 = b;
   while (v3 != 0)
     {
        register long q , t1, t2 ,t3;
        q = u3 / v3;
        t1 = u1 - q * v1; t2 = u2 - q * v2; t3 = u3 - q * v3;
        u1 = v1; u2 = v2; u3 = v3; v1 = t1; v2 = t2; v3 = t3;
     }
   u = u1; 
   v = u2;
   return d=u3;
//    return u3;
} 

int32& ZpzDom<Std32>::invext
  ( int32& u, const int32 a, const int32 b ) const
{
    register long u3;
    register long v1,v3;
    u = 1; u3 = a;
    v1 = 0; v3 = b;
    while (v3 != 0)
    {
	register long q, t1, t3;
	q = u3 / v3;
	t1 = u - q * v1; t3 = u3 - q * v3;
	u = v1; u3 = v3; v1 = t1; v3 = t3;
    }
    return (u3<0?u=-u:u);
}  


void ZpzDom<Std32>::Init() 
{
}

void ZpzDom<Std32>::End()
{
}

