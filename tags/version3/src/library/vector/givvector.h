#ifndef _GIV_VECTOR_H_
#define _GIV_VECTOR_H_
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/vector/givvector.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givvector.h,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description:
// Description of dense vector over T with classic arithmetic operations
// over T (vector x vector, vector x T, scalar product, shift)
// Vector handle computation over sub part of continuous elements of
// a vector as well as stride.

#include "givaro/giverror.h"
#include "givaro/givarray0.h"
#include "givaro/givcategory.h"
#include "givaro/giviterator.h"
#include "givaro/givops.h"


// ==========================================================================
// --
// -- VectorDom<Domain,StorageTag>:
// --
/** VectorDom<Domain,StorageTag>
*/
// ==========================================================================
template<class Domain, class StorageTag> class VectorDom {};

template<class Domain> class VectorDom<Domain,Dense>;
template<class Domain> class VectorDom<Domain,Sparse>;

#include "givaro/givvectorsparse.h"
#include "givaro/givvectordense.h"
// -- should be included in specialized class #include "givaro/givvectops.inl"
#include "givaro/givvectsparseops.inl"
#include "givaro/givvectdenseops.inl"

#endif
