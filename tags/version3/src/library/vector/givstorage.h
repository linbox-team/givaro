#ifndef _GIV_VECTOR_STORAGE_H_
#define _GIV_VECTOR_STORAGE_H_
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/vector/givstorage.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givstorage.h,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description:

#include "givaro/givcategory.h"

// ==========================================================================
// --
// -- RetVectorStorage<T, StorageTag>:
// -- return the storage type Storage_t associated with the
// -- StorageTag
// ==========================================================================
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



#endif
