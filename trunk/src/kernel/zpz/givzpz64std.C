#include <givaro/givconfig.h>
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givzpz64std.C,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givzpz64std.C,v 1.3 2004-10-11 13:54:38 jgdumas Exp $
// ==========================================================================
// Description:

#include <iostream>
#include "givaro/givzpz64std.h"

  // Returns d, and u and v such that u a + v b = d;
// JGD 04.11.1999
// const int64 ZpzDom<Std64>::gcdext 
//   ( int64& u, int64& v, const int64 a, const int64 b )
int64& ZpzDom<Std64>::gcdext 
  ( int64& d,  int64& u, int64& v, const int64 a, const int64 b ) const
{
   register int64 u1,u2,u3;
   register int64 v1,v2,v3;
   u1 = 1; u2 = 0; u3 = a;
   v1 = 0; v2 = 1; v3 = b;
   while (v3 != 0)
     {
        register int64 q , t1, t2 ,t3;
        q = u3 / v3;
        t1 = u1 - q * v1; t2 = u2 - q * v2; t3 = u3 - q * v3;
        u1 = v1; u2 = v2; u3 = v3; v1 = t1; v2 = t2; v3 = t3;
     }
   u = u1; 
   v = u2;
   return d=u3;
//    return u3;
} 

int64& ZpzDom<Std64>::invext 
  ( int64& u, const int64 a, const int64 b ) const
{
   register int64 u3;
   register int64 v1,v3;
   u = 1; u3 = a;
   v1 = 0; v3 = b;
   while (v3 != 0)
     {
        register int64 q , t1, t3;
        q = u3 / v3;
        t1 = u - q * v1; t3 = u3 - q * v3;
        u = v1; u3 = v3; v1 = t1; v3 = t3;
     }
   return (u3<0?u=-u:u); 
} 

void ZpzDom<Std64>::Init() 
{
} 

void ZpzDom<Std64>::End()
{
}
