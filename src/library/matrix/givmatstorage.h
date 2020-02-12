// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/matrix/givmatstorage.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givmatstorage.h,v 1.3 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================
// Description:


#ifndef _GIV_MATRIX_STORAGE_H_
#define _GIV_MATRIX_STORAGE_H_

#include "givaro/givcategory.h"

namespace Givaro {
    // =========================================================================
    // --
    // -- RetMatrixStorage<T, StorageTag>:
    // -- return the storage type Storage_t associated with the
    // -- StorageTag
    // =========================================================================
    template<class T, class StorageTag >
    struct RetMatrixStorage {
        typedef T		Type_t;
        typedef Undefined 	Storage_t;
    };


    // =========================================================================
    // --
    // -- RetMatrix2Storage<StorageTag,ViewTag>:
    // -- return the storage type associated with a view of the storage
    // -- associated with StorageTag.
    // =========================================================================
    template<class StorageTag, class ViewTag>
    struct RetMatrix2Storage {
        typedef Undefined Storage_t;
        typedef Undefined ViewTag_t;
    };

} // Givaro

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
