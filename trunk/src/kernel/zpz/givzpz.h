// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givzpz.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software. 
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givzpz.h,v 1.5 2009-09-17 14:28:23 jgdumas Exp $
// ==========================================================================
// Description:
//   Family of arithmetics over Zpz
#ifndef _GIVARO_ZPZ_H_ 
#define _GIVARO_ZPZ_H_ 

// ==========================================================================
// -- 
// ==========================================================================
#include <givaro/givconfig.h>

template<class TAG> class ZpzDom;

// -- Tag for arithmetic:
struct Std16 { typedef int16 type;}; // -- standard arithmetic over 16bits representations.
struct Std32 {typedef int32 type;}; // -- standard arithmetic over 32bits representations.
struct Unsigned32 {typedef uint32  type;}; // -- standard arithmetic over 32bits representations.

struct Log16 { typedef int16 type;}; // -- log arithmetic over 16bits representations.

#include "givaro/givzpz16std.h"
#include "givaro/givzpz16table1.h"
#include "givaro/givzpz32std.h"
#include "givaro/givzpz32uns.h"


#ifndef __USE_Givaro_SIXTYFOUR__
#ifdef __USE_64_bits__
#define __USE_Givaro_SIXTYFOUR__ 1
#endif

#ifdef __USE_ISOC99
#define __USE_Givaro_SIXTYFOUR__ 1
#endif
#endif

#ifndef __DONOTUSE_Givaro_SIXTYFOUR__
#ifndef __USE_GMPPLUSPLUS_SIXTYFOUR__
#define __USE_GMPPLUSPLUS_SIXTYFOUR__ 1
#endif
struct Std64 { typedef int64 type;}; // -- standard arithmetic over 64bits representations.
#include "givaro/givzpz64std.h"
#endif

#include "givaro/givzpzInt.h"

#endif
