// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/vector/givstorage.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givstorage.h,v 1.2 2009-09-17 14:28:23 jgdumas Exp $
// ==========================================================================
// Description:
#ifndef _GIV_VECTOR_STORAGE_H_
#define _GIV_VECTOR_STORAGE_H_

#include "givaro/givcategory.h"
namespace Givaro {

    // =========================================================================
    // --
    // -- RetVectorStorage<T, StorageTag>:
    // -- return the storage type Storage_t associated with the
    // -- StorageTag
    // =========================================================================
    template<class T, class StorageTag >
    struct RetVectorStorage {
        typedef T		Type_t;
        typedef Undefined 	Storage_t;
    };


    // ==========================================================================
    // --
    // -- RetVector2Storage<StorageTag,ViewTag>:
    // -- return the storage type associated with a view of the storage
    // -- associated with StorageTag.
    // ==========================================================================
    template<class StorageTag, class ViewTag>
    struct RetVector2Storage {
        typedef Undefined Storage_t;
        typedef Undefined ViewTag_t;
    };

} // Givaro

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
