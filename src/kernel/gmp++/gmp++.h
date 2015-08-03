// ========================================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Givaro version of gmp++.h
// Time-stamp: <03 Aug 15 11:12:00 Jean-Guillaume.Dumas@imag.fr>
// ========================================================================
#ifndef __GIVARO_GMPplusplus_H
#define __GIVARO_GMPplusplus_H
#include <string.h>
#include <limits.h>
#include <climits> // required by gcc 4.3
#include <givaro-config.h>


#ifdef __GIVARO_INLINE_ALL
#define giv_all_inlined inline
#else
#define giv_all_inlined
#endif

#ifndef __GIVARO__DONOTUSE_longlong__
#define __USE_64_bits__
#endif
#ifndef __DONOTUSE_Givaro_SIXTYFOUR__
#define __USE_64_bits__
#endif

#if !defined(GMP_VERSION_3) && !defined(GMP_NO_CXX) && !defined(__GIVARO_GMP_VERSION_3) && !defined(__GIVARO_GMP_NO_CXX)
// gmpxx.h defines __GMP_PLUSPLUS__
#include <gmpxx.h>
#endif

#ifdef __GIVARO_GMP_VERSION_3
extern "C" {
#endif

#include <gmp.h>

#ifdef __GIVARO_GMP_VERSION_3
}
#endif

#ifdef NDEBUG
#ifdef DEBUG
#error "NDEBUG and DEBUG both defined"
#endif
#endif

#ifndef _GIVARO_ISPRIMETESTS_
#define _GIVARO_ISPRIMETESTS_ 5
#endif

#include <gmp++/gmp++_int.h>
#include <gmp++/gmp++_rat.h>

#endif // __GIVARO_GMPplusplus_H
