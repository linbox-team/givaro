// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/matrix/givmatsparse.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givmatsparse.h,v 1.4 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================
// Description:
// of matrix by blocks.

#ifndef _GIV_MATRIX_SPARSE_H_
#define _GIV_MATRIX_SPARSE_H_

#include "givaro/givmatrix.h"
#include "givaro/givvector.h"
#include "givaro/givmatstoragesparse.h"

namespace Givaro {

    // --
    // -- Matrix class: dense matrix
    // --
    template <class Domain>
    class MatrixDom<Domain, Sparse>
    {
        Domain _domain;

    public:
        typedef 	   Domain 					Domain_t;
        typedef typename Domain::Element		Type_t;
        typedef 	   size_t					Indice_t;
        typedef 	   Dense 					StorageTag_t;
        typedef typename RetMatrixStorage<Type_t,Sparse>::Storage_t 	Storage_t;

        // -- Representation of Element of the domain
        typedef 	   Storage_t				Element;

        // -- Self_t
        typedef          MatrixDom<Domain, Sparse>			Self_t;

        //-- Dstor:
        ~MatrixDom() {}

        //-- Default cstor:
        MatrixDom() : _domain() {}

        //-- cstor:
        MatrixDom(const Domain& D) : _domain(D) {}

        //-- Cstor of recopy: compiler's generated
        MatrixDom(const Self_t& M )
        : _domain(M._domain) {}

        //-- init a new object, memory allocation
        void init(Element& r, Indice_t nr, Indice_t nc) const
        { r.allocate(nr, nc); }
        void init(Element& r) const
        { r.allocate(0,0); }
        void init(Element& A, const Element& B) const
        { A.copy(B); }

        //-- Assignment operator: physical copy
        void assign (Element& r, const Element& a) const
        { r.copy(a); }

        // -- Comparaizon
        int areEqual ( const Element& P, const Element& Q) const;
        int areNEqual( const Element& P, const Element& Q) const;
        int isZero  ( const Element& P ) const;

        //-- Dimension of the matrix space
        Indice_t nrow(const Element& A) const { return A._nrow; }
        Indice_t ncol(const Element& A) const { return A._ncol; }
        Domain_t subdomain() const { return _domain; }

        // -- arithmetic operator: operands could be aliased
        void mulin ( Element& res, const Type_t& u ) const;
        void mul   ( Element& res, const Type_t& u, const Element& v ) const;
        void mul   ( Element& res, const Element& u, const Type_t& v ) const;

        // VD is the vector domain for res and u
        void mul      ( typename VectorDom<Domain,Dense>::Element& res,
                        const Element& M,
                        const VectorDom<Domain,Dense>& VD,
                        const typename VectorDom<Domain,Dense>::Element& u ) const;
        void multrans ( typename VectorDom<Domain,Dense>::Element& res,
                        const Element& M,
                        const VectorDom<Domain,Dense>& VS,
                        const typename VectorDom<Domain,Dense>::Element& u ) const;


        void negin ( Element& P ) const
        {
            size_t sz = P._data.size();
            for(size_t i=0; i<sz; ++i) _domain.negin(P._data[i]);
        }

        void neg   ( Element& res, const Element& u ) const;

        // -- map of a unary operator, with operator()( Type_t& res)
        template<class OP>
        void map ( Element& res, OP& op ) const;

        // -- map of a unary operator, with operator()( Type_t& res, const Type_t& val)
        template<class OP>
        void map ( Element& res, OP& op, const Element& u ) const;

        // -- IO
        std::istream& read ( std::istream& s );
        std::ostream& write( std::ostream& s ) const;
        std::istream& read ( std::istream& s, Element& r ) const;
        std::ostream& write( std::ostream& s, const Element& r ) const;

        // -- Compression method to compact a dense matrix to a sparse
        // template<class StorageTag>,
        void compact( Element& Ms,
                      const MatrixDom<Domain, Dense>& MD,
                      const typename MatrixDom<Domain, Dense>::Element& Md);
    };

} // Givaro

//#include "givaro/givmatsparseops.inl"

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
