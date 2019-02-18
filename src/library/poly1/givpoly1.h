// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givpoly1.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givpoly1.h,v 1.3 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================
/** @file givpoly1.h
 * @ingroup poly1
 * @brief NO DOC
 */
#ifndef __GIVARO_poly1_H
#define __GIVARO_poly1_H

#include "givaro/giverror.h"
#include "givaro/givarray0.h"
#include "givaro/givcategory.h"
#include "givaro/giviterator.h"
#include "givaro/givops.h"

namespace Givaro {

    // ==========================================================================
    // --
    // -- Poly1Dom<Domain,StorageTag>:
    // --
    // ==========================================================================
    template<class Domain, class StorageTag=Dense> class Poly1Dom;

    template<class Domain> class Poly1Dom<Domain,Dense>;
    template<class Domain> class Poly1Dom<Domain,Sparse>;

} // Givaro

#include "givaro/givpoly1dense.h"
// #include "givaro/givvectorsparse.h"
// -- should be included in specialized class #include "givaro/givvectops.inl"
#include "givaro/givpoly1denseops.inl"
// #include "givaro/givvectsparseops.inl"

#endif // __GIVARO_poly1_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
