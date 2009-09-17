// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givpoly1.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software. 
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givpoly1.h,v 1.2 2009-09-17 14:28:23 jgdumas Exp $
// ==========================================================================
// Description:
#ifndef _GIV_POLY1_H_
#define _GIV_POLY1_H_

#include "givaro/giverror.h"
#include "givaro/givarray0.h"
#include "givaro/givcategory.h"
#include "givaro/giviterator.h"
#include "givaro/givops.h"


// ==========================================================================
// --
// -- Poly1Dom<Domain,StorageTag>:
// --
// ==========================================================================
template<class Domain, class StorageTag> class Poly1Dom;

template<class Domain> class Poly1Dom<Domain,Dense>;
template<class Domain> class Poly1Dom<Domain,Sparse>;

#include "givaro/givpoly1dense.h"
// #include "givaro/givvectorsparse.h"
// -- should be included in specialized class #include "givaro/givvectops.inl"
#include "givaro/givpoly1denseops.inl"
// #include "givaro/givvectsparseops.inl"

#endif
