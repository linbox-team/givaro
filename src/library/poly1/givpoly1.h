#ifndef _GIV_POLY1_H_
#define _GIV_POLY1_H_
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givpoly1.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givpoly1.h,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description:

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
