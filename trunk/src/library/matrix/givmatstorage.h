// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/matrix/givmatstorage.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software. 
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givmatstorage.h,v 1.2 2009-09-17 14:28:23 jgdumas Exp $
// ==========================================================================
// Description:
#ifndef _GIV_MATRIX_STORAGE_H_
#define _GIV_MATRIX_STORAGE_H_

#include "givaro/givcategory.h"

// ==========================================================================
// --
// -- RetMatrixStorage<T, StorageTag>:
// -- return the storage type Storage_t associated with the
// -- StorageTag
// ==========================================================================
template<class T, class StorageTag > 
struct RetMatrixStorage {
  typedef T		Type_t;
  typedef Undefined 	Storage_t;
};


// ==========================================================================
// --
// -- RetMatrix2Storage<StorageTag,ViewTag>:
// -- return the storage type associated with a view of the storage
// -- associated with StorageTag.
// ==========================================================================
template<class StorageTag, class ViewTag> 
struct RetMatrix2Storage {
  typedef Undefined Storage_t;
  typedef Undefined ViewTag_t; 
};



#endif
