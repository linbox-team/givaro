#ifndef _GIVARO_POINTER_H_
#define _GIVARO_POINTER_H_
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/memory/givpointer.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Author: T. Gautier
// $Id: givpointer.h,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description:
// - auto ptr management

#include "givaro/givaromm.h"


// ==================================================================== //

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

#endif
