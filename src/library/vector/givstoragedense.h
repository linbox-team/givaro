// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/vector/givstoragedense.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givstoragedense.h,v 1.2 2009-09-17 14:28:23 jgdumas Exp $
// ==========================================================================
// Description:
#ifndef _GIV_VECTOR_STORAGE_DENSE_H_
#define _GIV_VECTOR_STORAGE_DENSE_H_

#include "givaro/givstorage.h"

namespace Givaro {
    // =========================================================================
    // --
    // -- Specialization for dense representation, using Array0 of givaro
    // --
    // =========================================================================
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

} // Givaro
#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
