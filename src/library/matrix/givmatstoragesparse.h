// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/matrix/givmatstoragesparse.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givmatstoragesparse.h,v 1.3 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================
// Description:


#ifndef _GIV_MATRIX_STORAGE_SPARSE_H_
#define _GIV_MATRIX_STORAGE_SPARSE_H_

#include "givaro/givmatstorage.h"

namespace Givaro {
    // =========================================================================
    // --
    // -- Specialization for sparse representation, using Array0 of givaro
    // --
    // =========================================================================
    template<class T>
    struct RetMatrixStorage<T,Sparse> {
        typedef 		T			Type_t;
        typedef typename 	Array0<T>::Indice_t 	Indice_t;

        // --
        // -- Iterators of the container
        // --
        typedef typename Array0<T>::Iterator_t 	Iterator_t;
        typedef typename Array0<T>::constIterator_t 	constIterator_t;

        // --
        // --
        // --
        struct Storage_t {
            Indice_t   _nrow, _ncol;  // - dimension of the block matrix
            Array0<Indice_t>  _rows;  // - rows[i] of the first entry of
            // the i-th row in (_index,_data).
            // rows a +1 entry
            Array0<Indice_t>  _index; // - index of the first entry of
            Array0<T>         _data;  // - entries

            //-- Store dimension but don't allocate
            void allocate( size_t rsz, size_t csz )
            {
                _data.allocate(0);
                _index.allocate(0);
                _rows.allocate(rsz+1);
                _nrow = rsz; _ncol = csz;
            }
            //-- Store dimension but don't allocate
            void resize( size_t rsz, size_t csz )
            {
                _data.resize(0);
                _index.resize(0);
                _rows.resize(rsz+1);
                _nrow = rsz; _ncol = csz;
            }
            Storage_t& copy (const Storage_t& V)
            {
                _data.copy(V._data);
                _index.copy(V._index);
                _rows.copy(V._rows);
                _nrow = V._nrow; _ncol = V._ncol;
                return *this;
            }
            Storage_t& operator= (const Storage_t& V)
            {
                _data.copy(V._data);
                _index.copy(V._index);
                _rows.copy(V._rows);
                _nrow = V._nrow; _ncol = V._ncol;
                return *this;
            }
            Storage_t ( )
            : _nrow(0), _ncol(0), _rows(0), _index(0), _data( 0 ) {}
            Storage_t ( Indice_t nrow, Indice_t ncol )
            : _nrow(nrow), _ncol(ncol), _rows(nrow+1), _index(0), _data( 0 ) {}
            Storage_t ( const Storage_t& s )
            : _nrow(s._nrow), _ncol(s._ncol),
            _rows ( s._rows, givWithCopy()),
            _index( s._index, givWithCopy()),
            _data ( s._data, givWithCopy() ) {}
        };


    };

} // Givaro
#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
