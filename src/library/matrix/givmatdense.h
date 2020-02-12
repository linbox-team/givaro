// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/matrix/givmatdense.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givmatdense.h,v 1.4 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================
// Description:
// of matrix by blocks.
//

#ifndef _GIV_MATRIX_DENSE_H_
#define _GIV_MATRIX_DENSE_H_

#include "givaro/givmatrix.h"
#include "givaro/givmatstoragedense.h"


namespace Givaro {
    // --
    // -- Matrix class: dense matrix
    // --
    template <class Domain>
    class MatrixDom<Domain, Dense> {
        Domain 		  _domain;  		// -- domain of the entry
        VectorDom<Domain,Dense> _supportdomain;	// -- domain for some op.
    public:
        typedef 	   Domain 					Domain_t;
        typedef typename Domain::Element		Type_t;
        typedef 	   size_t					Indice_t;
        typedef 	   Dense 					StorageTag_t;
        typedef typename RetMatrixStorage<Type_t,Dense>::Storage_t    Storage_t;

        // -- Representation of Element of the domain
        typedef 	   Storage_t				Element;

        // -- Self_t
        typedef        MatrixDom<Domain, Dense>	Self_t;

        //-- Dstor:
        ~MatrixDom() {}

        //-- Default cstor:
        MatrixDom() : _domain(), _supportdomain() {}

        //-- Cstor of recopy: compiler's generated
        MatrixDom(const Self_t& M ) : _domain(M._domain), _supportdomain(M._supportdomain) {}

        // -- Cstor of a nr x nc dimensional matrix space
        MatrixDom(const Domain_t& D) : _domain(D), _supportdomain(D) {}

        //-- init a new object, memory allocation
        void init(Element& r, Indice_t nr, Indice_t nc) const
        { r.allocate(nr, nc); }

        void init(Element& r)
        { r.allocate(0); }

        //-- access operators:
        Type_t& operator() (Element& r, Indice_t i, Indice_t j) const
        { return r(i,j); }
        const Type_t& operator() (const Element& r, Indice_t i, Indice_t j) const
        { return r(i,j); }

        //-- Assignment operator: physical copy
        void assign (Element& r, const Element& a)
        { r.copy(a); }

        // -- Comparaizon
        int areEqual ( const Element& P, const Element& Q) const;
        int areNEqual( const Element& P, const Element& Q) const;
        int isZero  ( const Element& P ) const;

        //-- Dimension of the matrix space
        Indice_t nrow(const Element& r) const { return r.nrow(); }
        Indice_t ncol(const Element& r) const { return r.ncol(); }
        Domain_t subdomain() const { return _domain; }

        // -- arithmetic operators: operands could be aliased
        void mulin ( Element& res, const Element& u ) const;
        void mul   ( Element& res, const Element& u, const Element& v ) const;
        void addin ( Element& res, const Element& u ) const;
        void add   ( Element& res, const Element& u, const Element& v ) const;
        void subin ( Element& res, const Element& u ) const;
        void sub   ( Element& res, const Element& u, const Element& v ) const;
        void negin ( Element& res ) const;
        void neg   ( Element& res, const Element& u ) const;

        // --- Mul Vect:
        void mul   ( typename VectorDom<Domain,Dense>::Element& res, const Element& M,
                     const VectorDom<Domain,Dense>& VD,
                     const typename VectorDom<Domain,Dense>::Element& u ) const;
        void multrans ( typename VectorDom<Domain,Dense>::Element& res, const Element& M,
                        const VectorDom<Domain,Dense>& VS,
                        const typename VectorDom<Domain,Dense>::Element& u ) const;

        // -- axpy operations K-Space:
        // r <- a*x+y
        void axpy  ( Element& res, const Type_t& a, const Element& x, const Element& y )const;
        // r <- r+a*x
        void axpyin( Element& res, const Type_t& a, const Element& x ) const;
        // r <- y-a*x
        void axmy  ( Element& res, const Type_t& a, const Element& x, const Element& y ) const;
        // r <- r-a*x
        void axmyin( Element& res, const Type_t& a, const Type_t& x ) const;

        // a*A*X + bY
        void axpy  ( Element& res, const Type_t& a, const Element& A, const Element& X,
                     const Type_t& b, const Element& Y ) const;
        // A*X + Y
        void axpy  ( Element& res, const Element& A, const Element& X, const Element& Y ) const;

        // -- Element wise operation:
        void mulin ( Element& res, const Type_t& u ) const;
        void mul   ( Element& res, const Type_t& u, const Element& v ) const;
        void mul   ( Element& res, const Element& u, const Type_t& v ) const;

        // -- addition with a scalar: addition with I*val
        void add   ( Element& res, const Element& u, const Type_t& val ) const;
        void add   ( Element& res, const Type_t& val, const Element& v ) const;

        // -- substraction with a scalar: substraction with I*val
        void sub   ( Element& res, const Element& u, const Type_t& val ) const;
        void sub   ( Element& res, const Type_t& val, const Element& v ) const;

        // -- map of a inplace unary operator, with operator()( Type_t& res)
        template<class OP>
        void map ( Element& res, OP& op ) const;

        // -- map of a unary operator, with operator()( Type_t& res, const Type_t& val)
        template<class OP>
        void map ( Element& res, OP& op, const Element& u ) const;

        // -- with operator()( Type_t& res, const Type_t& v1, const Type_t& v2)
        template<class OP>
        void map ( Element&, OP&, const Element&, const Element&) const;

        // -- IO
        std::istream& read ( std::istream& s );
        std::ostream& write( std::ostream& s ) const;
        std::istream& read ( std::istream& s, Element& r ) const;
        std::ostream& write( std::ostream& s, const Element& r ) const;
    };

} // Givaro

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
