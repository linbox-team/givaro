#ifndef _GIV_PRIMES16_H_
#define _GIV_PRIMES16_H_
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givprimes16.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givprimes16.h,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description:
// - set of primes less than 2^16

#include <stddef.h>
#include "givconfig.h"

  // ---------------------------------------------  class Primes16
class Primes16  {
public:

  // -- Returns the number of primes of this ctxt
static size_t count() { return _size; } 

static int16 ith(int i) { 
   GIVARO_ASSERT( (i>=0)&&(i<(int)_size), "[Primes16::ith] index out of bounds");
   return _primes[i]; 
} 

private:
static const size_t _size;
static const int16  _primes[];
};


#endif
