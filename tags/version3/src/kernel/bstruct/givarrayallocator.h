#ifndef _ARRAY_ALLOCATOR_H_
#define _ARRAY_ALLOCATOR_H_
// ========================================================================== 
// $Source: /var/lib/cvs/Givaro/src/kernel/bstruct/givarrayallocator.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Author: T. Gautier
// $Id: givarrayallocator.h,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ========================================================================== 
// Description:


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
