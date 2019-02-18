// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/field/montgomery.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: A. Breust
// $Id: givzpz.h,v 1.8 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================

/*!@file montgomery.h
 * @ingroup zpz
 * @brief   Family of arithmetics over Zpz (\f$\mathbf{Z}/p\mathbf{Z}\f$).
 */

#ifndef __GIVARO_montgomery_H
#define __GIVARO_montgomery_H

// ==========================================================================
// --
// ==========================================================================

namespace Givaro {
    template<class TAG>
    class Montgomery;
}

#include "givaro/montgomery-int32.h"
#include "givaro/montgomery-ruint.h"

#endif // __GIVARO_montgomery_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
