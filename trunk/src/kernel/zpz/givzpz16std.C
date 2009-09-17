// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givzpz16std.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software. 
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givzpz16std.C,v 1.2 2009-09-17 14:28:23 jgdumas Exp $
// ==========================================================================
// Description:

#include <iostream>
#include "givaro/givzpz16std.h"

  // Returns d, and u and v such that u a + v b = d;
// JGD 04.11.1999
// const int32 ZpzDom<Std16>::gcdext 
//   ( int32& u, int32& v, const int32 a, const int32 b )
int32& ZpzDom<Std16>::gcdext 
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
} 

int32& ZpzDom<Std16>::invext 
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

void ZpzDom<Std16>::Init() 
{
} 

void ZpzDom<Std16>::End()
{
}

