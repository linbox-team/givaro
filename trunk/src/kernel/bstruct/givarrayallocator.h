// ========================================================================== 
// $Source: /var/lib/cvs/Givaro/src/kernel/bstruct/givarrayallocator.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software. 
// see the COPYRIGHT file for more details.
// Author: T. Gautier
// $Id: givarrayallocator.h,v 1.2 2009-09-17 14:28:22 jgdumas Exp $
// ========================================================================== 
#ifndef _ARRAY_ALLOCATOR_H_
#define _ARRAY_ALLOCATOR_H_

// -- ArrayAllocator: class for allocation of arrays. Should have
// - allocate(size_n)
// - reallocate(size_n)
// - destroy
template<class T, class Tag>
class ArrayAllocatort { };


// -- Specialization: for Array0Tag
class Array0Tag {};
template<class T, Array0Tag>
class ArrayAllocatort : public Array0<T> {};

#endif
