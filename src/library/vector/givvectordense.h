// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/vector/givvectordense.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givvectordense.h,v 1.3 2009-09-17 14:28:23 jgdumas Exp $
// ==========================================================================
// Description:
// Domain of dense vector over K with classic arithmetic operations
// over T (vector x vector, vector x T, scalar product, shift).
// A Element of this domain is a vector of K^n for any n.
// Vector handle computation over sub part of continuous Elements of
// a vector as well as stride.
#ifndef _GIV_VECTOR_DENSE_H_
#define _GIV_VECTOR_DENSE_H_

#include "givaro/givvector.h"
#include "givaro/givvectorsparse.h"
#include "givaro/givstoragedense.h"
#include "givaro/givelem.h"
namespace Givaro {

    template<class Domain>
    class VectorDom<Domain, Dense> {
    public:
        Domain 	_domain; // domain of the entries
    public:
        // -- Exported types
        typedef 	   Domain 					Domain_t;
        typedef typename Domain::Element 		Type_t;
        typedef 	   size_t					Indice_t;
        typedef 	   Dense	 				StorageTag_t;
        typedef typename RetVectorStorage<Type_t,Dense>::Storage_t 	Storage_t;

        // -- The representation of a dense vector
        typedef 	   Storage_t					Element;
        typedef typename VectorDom<Domain,Sparse>::Element SparseVector;
        // -- Self_t
        typedef 	   VectorDom<Domain, Dense> 			Self_t;

        // -- Dstor
        ~VectorDom() {}

        // -- Cstor of a new vector space of dimension s
        VectorDom<Domain,Dense>() : _domain() {}

        VectorDom<Domain,Dense>( const Domain& dom ) : _domain(dom) {}

        // -- Cstor of recopy
        VectorDom<Domain,Dense>( const Self_t& V ) : _domain(V._domain) {}

        bool operator==( const VectorDom<Domain,Dense>& BC) const
        { return _domain == BC._domain;}
        bool operator!=( const VectorDom<Domain,Dense>& BC) const
        { return _domain != BC._domain;}


        // -- init :
        void init( Element& v, size_t dim =0 ) const
        {
            v.resize(dim);
        }
        // --
        void init( Element& v, const Element& u ) const
        {
            v.copy(v);
        }

        // -- assignment operator: from a vector of the same vect space
        void assign ( Element& r, const Element& v) const
        {
            r.copy(v);
        }
        // -- assignment operator: from value * [1,.....1]
        void assign ( Element& r, size_t dim, const Type_t& val) const
        {
            r.resize( dim );
            for (size_t i=0; i<dim; ++i) r[i] = val;
        }

        // -- Comparizon
        bool areEqual (const Element& P, const Element& Q) const
        {
            size_t sP = P.size(), sQ = Q.size();
            if (sP != sQ) return false;
            for (size_t i=0; i<sP; ++i)
                if (!_domain.areEqual(P[i], Q[i])) return false;
            return true;
        }
        bool areNEqual(const Element& P, const Element& Q) const
        {
            return !areEqual(P,Q);
        }
        bool isZero(const Element& P) const
        {
            size_t sP =P.size();
            if (sP ==0) return true;
            for (size_t i=0; i<sP; ++i)
                if (!_domain.isZero(P[i])) return false;
            return true;
        }

        // -- Return the domain of the entries
        const Domain& subdomain() const { return _domain; }

        // -- return the dimension of the vector space
        size_t dim( const Element& v ) const { return v.size(); }

        // -- arithmetic operator: operands could be aliased
        void addin ( Element& res, const Element& u ) const;
        void add ( Element& res, const Element& u, const Element& v ) const;

        void addin( Element& res, const SparseVector& v ) const;
        void add ( Element& res, const Element& u, const SparseVector& v ) const;
        void add ( Element& res, const SparseVector& u, const Element& v ) const;

        void subin( Element& res, const Element& u ) const;
        void sub  ( Element& res, const Element& u, const Element& v ) const;

        void subin( Element& res, const SparseVector& v ) const;
        void sub ( Element& res, const Element& u, const SparseVector& v ) const;
        void sub ( Element& res, const SparseVector& u, const Element& v ) const;

        void negin ( Element& res ) const;
        void neg ( Element& res, const Element& u ) const;

        // - axpy like operations:
        // r <- a*x+y
        void axpy  ( Element& res, const Type_t& a, const Element& x, const Element& y )const;
        // r <- r+a*x
        void axpyin( Element& res, const Type_t& a, const Element& x ) const;
        // r <- y-a*x
        void axmy  ( Element& res, const Type_t& a, const Element& x, const Element& y ) const;
        // r <- r-a*x
        void axmyin( Element& res, const Type_t& a, const Element& x ) const;


        // Vector (+/-/*) Value ==  Element wise operation
        void mulin( Element& res, const Type_t& u ) const;
        void mul  ( Element& res, const Element& u, const Type_t& val ) const;
        void mul  ( Element& res, const Type_t& val, const Element& v ) const;

        void add ( Element& res, const Element& u, const Type_t& val ) const;
        void add ( Element& res, const Type_t& val, const Element& v ) const;

        void sub  ( Element& res, const Element& u, const Type_t& val ) const;
        void sub  ( Element& res, const Type_t& val, const Element& v ) const;


        // - axpy like operations, Element wise:
        // r <- a*x+y
        void axpy  ( Element& res, const Element& a, const Element& x, const Element& y ) const;
        // r <- r+a*x
        void axpyin( Element& res, const Element& a, const Element& x ) const;
        // r <- y-a*x
        void axmy  ( Element& res, const Element& a, const Element& x, const Element& y ) const;
        // r <- r-a*x
        void axmyin( Element& res, const Element& a, const Element& x ) const;


        // -- dot product: operands could be aliased
        void dot ( Type_t& res, const Element& u, const Element& v ) const;

        // -- map of a unary operator, with operator()( Type_t& res )
        // res and u could be aliases if OP permits it
        template<class UNOP>
        void map ( Element& res, UNOP& op ) const;

        // -- map of a unary operator, with operator()( Type_t& res, const Type_t& val)
        // res and u could be aliases if OP permits it
        template<class UNOP>
        void map ( Element& res, UNOP& op, const Element& u ) const;

        // -- map of a binary operator, with :
        // -- operator()( Type_t& res, const Type_t&, const Type_t& )
        template<class BINOP>
        void map ( Element&, const BINOP&, const Element&, const Element&) const;

        // -- IO
        std::istream& read ( std::istream& s );
        std::ostream& write( std::ostream& s ) const;
        std::istream& read ( std::istream& s, Element& r ) const;
        std::ostream& write( std::ostream& s, const Element& r ) const;

        // -- Iteration over a dense vector:
        typedef typename
        RetVectorStorage<Type_t,Dense>::Iterator_t          Iterator_t;
        typedef typename
        RetVectorStorage<Type_t,Dense>::constIterator_t     constIterator_t;
        typedef typename
        RetVectorStorage<Type_t,Dense>::IndiceIterator_t    IndiceIterator_t;

    };

} // Givaro

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
