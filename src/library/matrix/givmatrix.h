// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/matrix/givmatrix.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givmatrix.h,v 1.3 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================
// Description:
#ifndef _GIV_MATRIX_H_
#define _GIV_MATRIX_H_

#include "givaro/giverror.h"
#include "givaro/givarray0.h"
#include "givaro/givelem.h"

namespace Givaro {
    // =========================================================================
    // --
    // -- MatrixDom<T, StorageTag>:
    // -- current available tags are:
    // --  * Dense: dynamic n x m matrix
    // --  * Sparse: dynamic n x m matrix with few non null entries
    // --  * FixedBlock<n,m>: static n x m dimensional matrix
    // =========================================================================

    template <class T, class StorageTag> class MatrixDom { };

} // Givaro

#include "givaro/givvector.h"
#include "givaro/givmatstoragedense.h"
#include "givaro/givmatstoragesparse.h"
#include "givaro/givmatdense.h"
#include "givaro/givmatsparse.h"

#include "givaro/givmatdenseops.inl"
#include "givaro/givmatsparseops.inl"

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
