#ifndef _GIV_MATRIX_H_
#define _GIV_MATRIX_H_
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/matrix/givmatrix.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givmatrix.h,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description:

#include "givaro/giverror.h"
#include "givaro/givarray0.h"
#include "givaro/givelem.h"


// ==========================================================================
// --
// -- MatrixDom<T, StorageTag>:
// -- current available tags are:
// --  * Dense: dynamic n x m matrix
// --  * Sparse: dynamic n x m matrix with few non null entries
// --  * FixedBlock<n,m>: static n x m dimensional matrix
// ==========================================================================

template <class T, class StorageTag> class MatrixDom { };

#include "givaro/givvector.h"
#include "givaro/givmatstoragedense.h"
#include "givaro/givmatstoragesparse.h"
#include "givaro/givmatdense.h"
#include "givaro/givmatsparse.h"

#include "givaro/givmatdenseops.inl"
#include "givaro/givmatsparseops.inl"

#endif
