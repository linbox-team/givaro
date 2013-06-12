// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givzpz.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givzpz.h,v 1.8 2011-02-02 16:23:56 bboyer Exp $
// ==========================================================================

/*!@file givzpz.h
 * @ingroup zpz
 * @brief   Family of arithmetics over Zpz (\f$\mathbf{Z}/p\mathbf{Z}\f$).
 */

#ifndef __GIVARO_zpz_H
#define __GIVARO_zpz_H

// ==========================================================================
// --
// ==========================================================================
#include <givaro/givconfig.h>

namespace Givaro {
	template<class TAG>
	class ZpzDom;
}

#include "givaro/givzpztypes.h"

#include <givaro/givzpzGen.h>

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

#endif // __GIVARO_zpz_H
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
