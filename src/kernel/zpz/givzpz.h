// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givzpz.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givzpz.h,v 1.7 2011-02-02 13:45:03 jgdumas Exp $
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

#include "givaro/givzpztypes.h"

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
#include "givaro/givzpz64std.h"
#endif

#include "givaro/givzpzInt.h"

#endif
