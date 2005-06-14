#ifndef _GIV_ARRAY0_H_
#define _GIV_ARRAY0_H_
// ========================================================================== 
// $Source: /var/lib/cvs/Givaro/src/kernel/bstruct/givarray0.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Author: T. Gautier
// $Id: givarray0.h,v 1.3 2005-06-14 14:53:14 pernet Exp $
// ========================================================================== 
// Description:
// Array of type T with reference mecanism.
#include <stddef.h> // size_t

#include "givaro/givaromm.h"
#include "givaro/givperf.h"
#include "givaro/giverror.h"

GIVARO_PERF_DEFCLASS(Array0,T)

template <class T>
class Array0 GIVARO_PERF_INEHERIT(Array0,T) {
  //-- Recopy cstor : private: don't use it
  Array0 (const Array0<T>& p);
  void build( size_t s, const T& t) ;
public :
  typedef int  		Indice_t;
  typedef T    		Type_t;
  typedef Array0<T> 	Self_t;
  typedef Type_t	*Iterator_t;
  typedef const Type_t	*constIterator_t;

      // STL compliance
  typedef Type_t value_type;
  typedef Type_t *iterator;
  typedef const Type_t *const_iterator;


  //-- Default cstor : ctsor of s size array
  Array0 (size_t  s = 0);
  Array0 (size_t  s, const T& t);

  //-- Recopy cstor : logical copy
  Array0 (const Self_t& p, givNoCopy);
  //-- Recopy cstor : physical copy
  Array0 (const Self_t& p, givWithCopy);

  //-- Destructor
  ~Array0 ();

  //-- Destroy of the array
  void destroy ();

  //-- Allocation of an array of s Elements: if refcount>1
  // then it is always a creation of new array
  void allocate (size_t s);

  //-- Reallocation of an array of s Elements: if refcount>1
  // then it is always a creation of new array + recopy 
  void reallocate (size_t s);

  //-- Physical copy operator: reallocate dest of the same size
  // as src (if necessary) and apply GivaroCopyItem<Array<T>,T> on each Element.
  // This class can be specialized. Return dest (i.e, *this).
  Self_t& copy(const Self_t& src);

  //-- Logical recopy operator: make an alias to src. Return dest.
  Self_t& logcopy(const Self_t& src);

  //-- Return the occuped size of the array
  size_t size() const { return _size; }

  //-- Return the physical size of the array (capacity)
  size_t phsize() const;

  //-- Return the base ptr to the array 
  Type_t* baseptr();
  Type_t* const baseptr() const;

  //-- Access to the ith Element:
  const T& operator[] (Indice_t i)  const; //  { return _d[i]; }
  T& operator[] (Indice_t i); //  { return _d[i]; } ;
  void write(Indice_t i, const Type_t& val);
  void read (Indice_t i, Type_t& val) const;


  Iterator_t begin();
  Iterator_t end();
  constIterator_t begin() const;
  constIterator_t end() const;

protected :  //--------------------- protected Internal representation
  int* _cnt;     // reference counter on _d
  size_t _size;  // actual size of the array. If ==0 then _psz=_d=_cnt=0
  size_t _psz;   // physical size of the array
  T*  _d;        // ptr to the memory
private:
  //-- assignement operator cannot be herited. 
  Self_t& operator= (const Self_t& p);
};


#include "givaro/givarray0.inl"

#endif
