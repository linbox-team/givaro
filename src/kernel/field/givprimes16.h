// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givprimes16.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// Time-stamp: <02 Oct 07 16:43:17 Jean-Guillaume.Dumas@imag.fr>
// ==========================================================================

/*! @file givprimes16.h
 * @ingroup zpz
 * @brief  set of primes less than 2^16
 */

#ifndef __GIVARO_primes16_H
#define __GIVARO_primes16_H

#include <stddef.h>
namespace Givaro {
  //! class Primes16
class Primes16  {
public:

  // -- Returns the number of primes of this ctxt
static size_t count() { return _size; }

static size_t ith(size_t i) {
//    JGD 02.10.2007 : dependent of system/givconfig.h
//    GIVARO_ASSERT( (i>=0)&&(i<(int)_size), "[Primes16::ith] index out of bounds");
    return _primes[i];
}

private:
static const size_t _size;
static const size_t  _primes[];
};
} // namespace Givaro

#endif // __GIVARO_primes16_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
