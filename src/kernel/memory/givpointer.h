// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/memory/givpointer.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Author: T. Gautier
// $Id: givpointer.h,v 1.3 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================
/*! @file givpointer.h
 * @ingroup memory
 * @brief  auto ptr management
 */
#ifndef __GIVARO_pointer_H
#define __GIVARO_pointer_H

#include "givaro/givaromm.h"

namespace Givaro {

    // ==================================================================== //

    //! Refcount Pointer
    template<class T>
    class RefCountPtr {
        T*   _data;
        mutable int* _count;
    public:
        explicit RefCountPtr ( T* data )
        : _data( data ), _count(0)
        {
            _count = GivaroMM<int>::allocate(1);
            *_count = 1;
        }
        RefCountPtr ( const RefCountPtr<T>& ptr )
        : _data( ptr._data), _count(ptr._count)
        {
            if (_count !=0) *_count += 1;
        }
        ~RefCountPtr()
        {
            if (--*_count ==0) {
                delete data;
                GivaroMM<int>::desallocate(_count);
            }
        }

        RefCountPtr<T>& operator=( const RefCountPtr<T>& ptr )
        {
            if (--*_count ==0) {
                delete data;
                GivaroMM<int>::desallocate(_count);
            }
            _data = ptr._data; _count = ptr._count;
            if (_count !=0) *_count += 1;
        }

        T& operator* () const { return *_data; }
        T* operator-> () const { return _data; }
    };

} // namespace Givaro

#endif // __GIVARO_pointer_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
