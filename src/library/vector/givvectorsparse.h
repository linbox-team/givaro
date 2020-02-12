// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/vector/givvectorsparse.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givvectorsparse.h,v 1.3 2009-09-17 14:28:23 jgdumas Exp $
// ==========================================================================
// Description:
// Description of sparse vector over T with classic arithmetic operations
// over T (vector x vector, vector x T, scalar product).
// The storage for sparse vector are returned by
// - RetVectorStorage<T,Sparse>::Storage_t: that provide following
// properties:
// - iterator on the index domain iterates on increasing value
// -
#ifndef _GIV_VECTOR_SPARSE_H_
#define _GIV_VECTOR_SPARSE_H_

#include "givaro/givvector.h"
#include "givaro/givstoragesparse.h"
#include "givaro/givelem.h"
namespace Givaro {

    template<class Domain>
    class VectorDom<Domain, Sparse> {
        Domain _domain;	// domain of the entries
    public :
        // -- Exported types
        typedef typename Domain::Element		Type_t;
        typedef 	   Domain   				Domain_t;
        typedef 	   size_t					Indice_t;
        typedef 	   Sparse	 			 	StorageTag_t;
        typedef typename RetVectorStorage<Type_t,Sparse>::Storage_t 	Storage_t;

        // -- Representation of Element of VectorDom<D, Sparse>
        typedef 	   Storage_t				Element;

        // -- Self_t
        typedef 	   VectorDom<Domain, Sparse> 	 		Self_t;

        // -- Dstor
        ~VectorDom() {}

        // -- Cstor of a new vector of size s Elements
        VectorDom( const Domain& D = Domain() ) : _domain(D) {}

        // -- Cstor of recopy
        VectorDom(const Self_t& V) : _domain(V._domain) {}

        int operator==( const VectorDom<Domain,Sparse>& BC) const
        { return _domain == BC._domain;}
        int operator!=( const VectorDom<Domain,Sparse>& BC) const
        { return _domain != BC._domain;}

        // -- assignment operator: from a vector
        void init ( Element& r, size_t dim =0) const
        { r.allocate(dim,0); }

        // -- assignment operator: from a vector
        void assign (Element& r, const Element& v)
        {
            r.copy(v);
        }

        // -- Comparaizon
        int areEqual ( const Element& P, const Element& Q) const;
        int areNEqual( const Element& P, const Element& Q) const;
        int isZero  ( const Element& P ) const;

        // -- return the dimension of a vector
        size_t dim( const Element& u ) const { return u.size(); }
        const Domain& subdomain() const { return _domain; }

        // -- Arithmetic operations: base
        void add ( Element& res, const Element& op1, const Element& op2) const;
        void sub ( Element& res, const Element& op1, const Element& op2) const;

        // -- dot product: operands could be aliased
        void dot ( Type_t& res, const Element& u, const Element& v ) const;

        // -- Syntaxic sugar: (Value) op (Vector): Element wise ops.
        void addin( Element& res, const Element& u ) const;
        void add  ( Element& res, const Element& u, const Type_t& val ) const;
        void add  ( Element& res, const Type_t& val, const Element& v ) const;
        void subin( Element& res, const Element& u ) const;
        void sub  ( Element& res, const Element& u, const Type_t& val ) const;
        void sub  ( Element& res, const Type_t& val, const Element& v ) const;
        void negin( Element& res ) const;
        void neg  ( Element& res, const Element& u ) const;

        // -- Compression method to compact a dense vector
        void compact( Element& u, const VectorDom<Domain, Dense>& VDom,
                      const typename VectorDom<Domain, Dense>::Element& v ) const;

        // -- Compression method to compact a sparse vector
        void compact( Element& u, const VectorDom<Domain, Sparse>& VDom,
                      const typename VectorDom<Domain, Sparse>::Element& v ) const;

        template<class UNOP>
        void map( Element& r, const UNOP& op, const Element& u) const;

        template<class UNOP>
        void map( Element& r, UNOP& op, const Element& u) const;

        // -- IO: domain
        std::ostream& write( std::ostream& o ) const;
        std::istream& read ( std::istream& i );

        // -- IO: domain Element
        std::ostream& write( std::ostream& o, const Element& r ) const;
        std::istream& read ( std::istream& i, Element& r ) const;


        // -- Iteration over a sparse vector:
        typedef typename RetVectorStorage<Type_t,Sparse>::Iterator_t 		Iterator_t;
        typedef typename RetVectorStorage<Type_t,Sparse>::constIterator_t 	constIterator_t;
        typedef typename RetVectorStorage<Type_t,Sparse>::IndiceIterator_t 	IndiceIterator_t;

        Iterator_t	  begin_data( Element& U ) const { return U.begin_data(); }
        Iterator_t	  end_data  ( Element& U ) const { return U.end_data(); }
        constIterator_t begin_data( const Element& U ) const { return U.begin_data(); }
        constIterator_t end_data  ( const Element& U ) const { return U.end_data(); }
        IndiceIterator_t begin_indice( const Element& U ) const { return U.begin_indice(); }
        IndiceIterator_t end_indice  ( const Element& U ) const { return U.end_indice(); }
    };

} //Givaro
#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
