#ifndef _GIV_ARRAY_FIXED_H_
#define _GIV_ARRAY_FIXED_H_
// ========================================================================== 
// $Source: /var/lib/cvs/Givaro/src/kernel/bstruct/givarrayfixed.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Author: T. Gautier
// $Id: givarrayfixed.h,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ========================================================================== 
// Description:
// ArrayFixed of type T with fixed dimension
#include <stddef.h> // size_t

#include "givaro/givaromm.h"
#include "givaro/giverror.h"
#include "givaro/givperf.h"

template <class T, size_t SIZE>
class ArrayFixed GIVARO_PERF_INEHERIT(ArrayFixed,T) {
  T  _data[SIZE];        // _data
public :
  typedef int  			Indice_t;
  typedef T    			Type_t;
  typedef ArrayFixed<T,SIZE> 	Self_t;
  typedef Type_t		*Iterator_t;
  typedef const Type_t		*constIterator_t;

  //-- Default cstor: 
  ArrayFixed () {};

  //-- Recopy cstor : implicit

  //-- Destructor: implicit

  //-- Physical copy operator: reallocate dest of the same SIZE
  // as src (if necessary) and apply GivaroCopyItem<Array<T>,T> on each element.
  // This class can be specialized. Return dest (i.e, *this).
  Self_t& copy(const Self_t& src);

  //-- Return the occuped SIZE of the array
  size_t size() const { return SIZE; }

  //-- Return the physical size of the array (capacity)
  size_t phsize() const { return SIZE; }

  //-- Return the base ptr to the array 
  Type_t* baseptr() { return _data; }
  Type_t* const baseptr() const { return _data; }

  //-- Access to the ith element:
  const T& operator[] (Indice_t i)  const { 
    GIVARO_ASSERT((i >=0)&&(i<SIZE), "[Array<T>::[]]: index out of bounds.");
    return _data[i]; 
  }

  T& operator[] (Indice_t i) { 
    GIVARO_ASSERT((i >=0)&&(i<SIZE), "[Array<T>::[]]: index out of bounds.");
    return _data[i]; 
  }

  // -- Iterator
  Iterator_t begin() { return _data; }
  Iterator_t end() { return _data + SIZE; }
  constIterator_t begin() const { return _data; }
  constIterator_t end() const { return _data + SIZE; }

  template<class UNARYOP>
  void map( UNARYOP& opcode );

  template<class UNARYOP>
  void map( UNARYOP& opcode ) const;

protected :  //--------------------- protected Internal representation
private:
  //-- assignement operator cannot be herited. 
  Self_t& operator= (const Self_t& p) {};
};

// -- Map opcode on all elements less or requal that ith
// -- Terminal recursion, specialization
template<class T, class UNARYOP, size_t ith>
struct __giv_map_less_ith;

template<class T, class UNARYOP>
struct __giv_map_less_ith<T, UNARYOP, 0> {
  inline void operator()( register T* data, UNARYOP& opcode )
  { opcode(data[0]); }
};

template<class T, class UNARYOP, size_t ith>
struct __giv_map_less_ith<T, UNARYOP, ith> {
  inline void operator()( register T* data, UNARYOP& opcode )
  {  
    opcode(data[ith]); 
    __giv_map_less_ith<T,UNARYOP,ith-1>()(data, opcode); 
  }
};

template<class T, class UNARYOP, size_t ith>
struct __giv_map_less_ith_const;

template<class T, class UNARYOP>
struct __giv_map_less_ith_const<T,UNARYOP,0> {
  inline operator()( register const T* data, UNARYOP& opcode )
  { opcode(data[0]); }
};

template<class T, class UNARYOP, size_t ith>
struct __giv_map_less_ith_const<T,UNARYOP,ith> {
  inline void operator() ( register const T* data, UNARYOP& opcode )
  {  opcode(data[ith]); 
     __giv_map_less_ith<T,UNARYOP,ith-1>(data, opcode); 
  }
};


template<class T, size_t SIZE>
template<class UNARYOP>
void ArrayFixed<T,SIZE>::map( UNARYOP& opcode )
{ __giv_map_less_ith<T,UNARYOP,SIZE>()(_data, opcode); }

template<class T, size_t SIZE>
template<class UNARYOP>
void ArrayFixed<T,SIZE>::map( UNARYOP& opcode ) const
{ __giv_map_less_ith_const<T,UNARYOP,SIZE>()(_data, opcode); }



#endif
