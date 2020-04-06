// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/vector/givstoragesparse.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givstoragesparse.h,v 1.2 2009-09-17 14:28:23 jgdumas Exp $
// ==========================================================================
// Description:
#ifndef _GIV_VECTOR_STORAGE_SPARSE_H_
#define _GIV_VECTOR_STORAGE_SPARSE_H_

#include "givaro/givstorage.h"

namespace Givaro {
    // =========================================================================
    // --
    // -- Specialization for sparse representation, using Array0 over (int,T)
    // --
    // =========================================================================

    template<class T>
    struct RetVectorStorage<T,Sparse> {
        typedef T             Type_t;

        // --
        // -- Iterators
        // --
        typedef typename Array0<T>::Indice_t 	  Indice_t;
        typedef typename Array0<T>::Iterator_t  Iterator_t;
        typedef typename Array0<T>::constIterator_t constIterator_t;
        typedef typename Array0<Indice_t>::constIterator_t  IndiceIterator_t;

        // --
        // -- wrapper for Array<(I1,I2)> == (Array<I1>,Array<I2>)
        // --
        struct Storage_t {
            size_t 		_dim;
            Array0<Indice_t> 	_index;
            Array0<T>   	_data;

            size_t dim() const { return _dim; }
            size_t size() const { return _index.size(); }
            void allocate( size_t dim, size_t sz =0)
            { _dim = dim; _index.allocate(sz); _data.allocate(sz); }
            void resize( size_t dim, size_t sz =0)
            { _dim = dim; _index.resize(sz); _data.resize(sz); }
            Storage_t& copy (const Storage_t& V)
            {
                _index.copy(V._index);
                _data.copy(V._data);
                return *this;
            }
            Storage_t& operator= (const Storage_t& V)
            {
                _index.copy(V._index);
                _data.copy(V._data);
                return *this;
            }
            Storage_t ( const Array0<Indice_t>& i, Array0<T>& d  )
            : _index(i, givWithCopy() ), _data( d, givWithCopy() ) {}
            Storage_t ( const Storage_t& s )
            : _index(s._index, givWithCopy() ), _data( s._data, givWithCopy() ) {}
            Storage_t ( size_t sz =0 )
            : _index(sz), _data( sz ) {}

            typename Array0<T>::Iterator_t begin_data() { return _data.begin(); }
            typename Array0<T>::Iterator_t end_data()   { return _data.end(); }
            typename Array0<T>::constIterator_t begin_data() const { return _data.begin(); }
            typename Array0<T>::constIterator_t end_data() const   { return _data.end(); }
            typename Array0<Indice_t>::constIterator_t begin_indice() const { return _index.begin(); }
            typename Array0<Indice_t>::constIterator_t end_indice() const   { return _index.end(); }
        };

    };

} // Givaro

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
