// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/vector/givvector.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givvector.h,v 1.3 2009-09-17 14:28:23 jgdumas Exp $
// ==========================================================================
// Description:
// Description of dense vector over T with classic arithmetic operations
// over T (vector x vector, vector x T, scalar product, shift)
// Vector handle computation over sub part of continuous Elements of
// a vector as well as stride.
#ifndef _GIV_VECTOR_H_
#define _GIV_VECTOR_H_

#include "givaro/giverror.h"
#include "givaro/givarray0.h"
#include "givaro/givcategory.h"
#include "givaro/giviterator.h"
#include "givaro/givops.h"

namespace Givaro {
    // =========================================================================
    // --
    // -- VectorDom<Domain,StorageTag>:
    // --
    // =========================================================================
    template<class Domain, class StorageTag> class VectorDom {};

    template<class Domain> class VectorDom<Domain,Dense>;
    template<class Domain> class VectorDom<Domain,Sparse>;

} // Givaro

#include "givaro/givvectorsparse.h"
#include "givaro/givvectordense.h"
// -- should be included in specialized class #include "givaro/givvectops.inl"
#include "givaro/givvectsparseops.inl"
#include "givaro/givvectdenseops.inl"

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
