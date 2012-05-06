// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/bstruct/givarrayallocator.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Author: T. Gautier
// $Id: givarrayallocator.h,v 1.3 2011-02-02 16:23:55 bboyer Exp $
// ==========================================================================

/** @file givarrayallocator.h
 * @ingroup bstruct
 * @brief NO DOC
 */
#ifndef __GIVARO_array_allocator_H
#define __GIVARO_array_allocator_H

namespace Givaro {


/*! @brief ArrayAllocator: class for allocation of arrays.
 * Should have
 * - allocate(size_n)
 * - reallocate(size_n)
 * - destroy
 * .
 */
template<class T, class Tag>
class ArrayAllocatort { };


//! Array0Tag
class Array0Tag {};

//! Specialization: for Array0Tag
template<class T, Array0Tag>
class ArrayAllocatort : public Array0<T> {};

} // namespace Givaro


#endif // __GIVARO_array_allocator_H
