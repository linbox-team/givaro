// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/matrix/givmatstoragedense.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givmatstoragedense.h,v 1.3 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================
// Description:


#ifndef _GIV_MATRIX_STORAGE_DENSE_H_
#define _GIV_MATRIX_STORAGE_DENSE_H_

#include "givaro/givmatstorage.h"

namespace Givaro {
    // =========================================================================
    // --
    // -- Specialization for dense representation, using Array0 of givaro
    // --
    // =========================================================================
    template<class T>
    struct RetMatrixStorage<T,Dense> {
        typedef 	   T			Type_t;
        typedef typename Array0<T>::Indice_t 	Indice_t;

        // --
        // -- Iterators on the storage object: linearization
        // --
        typedef typename Array0<T>::Iterator_t 	Iterator_t;
        typedef typename Array0<T>::constIterator_t 	constIterator_t;

        // --
        // -- Storage: row storage organisation
        // --
        struct Storage_t : public Array0<T> {
            Indice_t _nrow;
            Indice_t _ncol;
            void allocate  ( Indice_t nrow, Indice_t ncol)
            { Array0<T>::allocate( nrow*ncol );
                _nrow = nrow; _ncol = ncol;
            }
            Type_t& operator() (Indice_t i, Indice_t j)
            { return Array0<T>::operator[]( i*_ncol+j ); }
            const Type_t& operator() (Indice_t i, Indice_t j) const
            { return Array0<T>::operator[]( i*_ncol+j ); }
            void resize( Indice_t nrow, Indice_t ncol)
            {
                Storage_t tmp; tmp.allocate(nrow, ncol);
                Indice_t i, mrow = (nrow < _nrow ? nrow : _nrow);
                Indice_t j, mcol = (ncol < _ncol ? ncol : _ncol);
                for (i=0; i<mrow; ++i)
                    for (j=0; j<mrow; ++j)
                        tmp(i,j) = (*this)(i,j);
                this->logcopy(tmp);
                _nrow = nrow; _ncol = ncol;
            };
            Indice_t nrow() const { return _nrow; }
            Indice_t ncol() const { return _ncol; }
        };
    };

} // givaro

#ifdef GIVARO_USE_SPECIALISED
#ifdef GIVARO_HAVE_LBLAS // -- specialization
// #include "givaro/givmatstoragedense.f.spe" // float
// #include "givaro/givmatstoragedense.d.spe" // double
#endif
#endif

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
