#ifndef _GIV_VECTOR_STORAGE_DENSE_H_
#define _GIV_VECTOR_STORAGE_DENSE_H_
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/vector/givstoragedense.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givstoragedense.h,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description:

#include "givaro/givstorage.h"


// ==========================================================================
// --
// -- Specialization for dense representation, using Array0 of givaro
// --
// ==========================================================================
template<class T> 
struct RetVectorStorage<T,Dense> {
  typedef T		Type_t; 

  // --
  // -- Iterators 
  // --
  typedef typename Array0<T>::Indice_t 		Indice_t;
  typedef typename Array0<T>::Iterator_t 	Iterator_t;
  typedef typename Array0<T>::constIterator_t 	constIterator_t;
  typedef 	   Array0<T> 			Storage_t;

  struct IndiceIterator_t {
  public:
    typedef IndiceIterator_t Self_t;

    // - cstor
    IndiceIterator_t(Indice_t a) : _curr(a) {};

    Indice_t operator*() const { return _curr; }
    Self_t&  operator++() { ++_curr; return *this;}
    Self_t&  operator--() { --_curr; return *this;}
    Self_t   operator++(int) { Self_t tmp = *this; ++_curr; return tmp; }
    Self_t   operator--(int) { Self_t tmp = *this; --_curr; ; return tmp; }
    Self_t&  operator+=(Indice_t n) { _curr += n; return *this; }
    Self_t&  operator-=(Indice_t n) { _curr -= n; return *this; }
    Indice_t operator[](Indice_t n) const { return _curr + n; }
    int operator==(const Self_t& n) const { return (_curr == n._curr); }
    int operator!=(const Self_t& n) const { return (_curr != n._curr); }

    Indice_t _curr;
  };
};


#endif
