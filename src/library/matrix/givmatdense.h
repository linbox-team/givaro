// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/matrix/givmatdense.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givmatdense.h,v 1.4 2011-02-02 16:23:56 bboyer Exp $
// ==========================================================================
// Description:
// of matrix by blocks.
//

#error "dead code"

#ifndef _GIV_MATRIX_DENSE_H_
#define _GIV_MATRIX_DENSE_H_

#include "givaro/givmatrix.h"
#include "givaro/givmatstoragedense.h"


// --
// -- Matrix class: dense matrix
// --
template <class Domain>
class MatrixDom<Domain, Dense> {
  Domain 		  _domain;  		// -- domain of the entry
  VectorDom<Domain,Dense> _supportdomain;	// -- domain for some op.
public:
  typedef 	   Domain 					Domain_t;
  typedef typename Domain::Rep					Type_t;
  typedef 	   int 						Indice_t;
  typedef 	   Dense 					StorageTag_t;
  typedef typename RetMatrixStorage<Type_t,Dense>::Storage_t    Storage_t;

  // -- Representation of Element of the domain
  typedef 	   Storage_t    				Rep;

  // -- Self_t
  typedef          MatrixDom<Domain, Dense>                     Self_t;

  //-- Dstor:
  ~MatrixDom() {}

  //-- Default cstor:
  MatrixDom() : _domain(), _supportdomain() {}

  //-- Cstor of recopy: compiler's generated
  MatrixDom(const Self_t& M ) : _domain(M._domain), _supportdomain(M._supportdomain) {}

  // -- Cstor of a nr x nc dimensional matrix space
  MatrixDom(const Domain_t& D) : _domain(D), _supportdomain(D) {}

  //-- init a new object, memory allocation
  void init(Rep& r, Indice_t nr, Indice_t nc) const
  { r.allocate(nr, nc); }

  void init(Rep& r)
  { r.allocate(0); }

  //-- access operators:
  Type_t& operator() (Rep& r, Indice_t i, Indice_t j) const
  { return r(i,j); }
  const Type_t& operator() (const Rep& r, Indice_t i, Indice_t j) const
  { return r(i,j); }

  //-- Assignment operator: physical copy
  void assign (Rep& r, const Rep& a)
  { r.copy(a); }

  // -- Comparaizon
  int areEqual ( const Rep& P, const Rep& Q) const;
  int areNEqual( const Rep& P, const Rep& Q) const;
  int iszero  ( const Rep& P ) const;

  //-- Dimension of the matrix space
  Indice_t nrow(const Rep& r) const { return r.nrow(); }
  Indice_t ncol(const Rep& r) const { return r.ncol(); }
  Domain_t subdomain() const { return _domain; }

  // -- arithmetic operators: operands could be aliased
  void mulin ( Rep& res, const Rep& u ) const;
  void mul   ( Rep& res, const Rep& u, const Rep& v ) const;
  void addin ( Rep& res, const Rep& u ) const;
  void add   ( Rep& res, const Rep& u, const Rep& v ) const;
  void subin ( Rep& res, const Rep& u ) const;
  void sub   ( Rep& res, const Rep& u, const Rep& v ) const;
  void negin ( Rep& res ) const;
  void neg   ( Rep& res, const Rep& u ) const;

  // --- Mul Vect:
  void mul   ( typename VectorDom<Domain,Dense>::Rep& res, const Rep& M,
               const VectorDom<Domain,Dense>& VD,
               const typename VectorDom<Domain,Dense>::Rep& u ) const;
  void multrans ( typename VectorDom<Domain,Dense>::Rep& res, const Rep& M,
                  const VectorDom<Domain,Dense>& VS,
                  const typename VectorDom<Domain,Dense>::Rep& u ) const;

  // -- axpy operations K-Space:
  // r <- a*x+y
  void axpy  ( Rep& res, const Type_t& a, const Rep& x, const Rep& y )const;
  // r <- r+a*x
  void axpyin( Rep& res, const Type_t& a, const Rep& x ) const;
  // r <- y-a*x
  void axmy  ( Rep& res, const Type_t& a, const Rep& x, const Rep& y ) const;
  // r <- r-a*x
  void axmyin( Rep& res, const Type_t& a, const Type_t& x ) const;

  // a*A*X + bY
  void axpy  ( Rep& res, const Type_t& a, const Rep& A, const Rep& X,
               const Type_t& b, const Rep& Y ) const;
  // A*X + Y
  void axpy  ( Rep& res, const Rep& A, const Rep& X, const Rep& Y ) const;

  // -- Element wise operation:
  void mulin ( Rep& res, const Type_t& u ) const;
  void mul   ( Rep& res, const Type_t& u, const Rep& v ) const;
  void mul   ( Rep& res, const Rep& u, const Type_t& v ) const;

  // -- addition with a scalar: addition with I*val
  void add   ( Rep& res, const Rep& u, const Type_t& val ) const;
  void add   ( Rep& res, const Type_t& val, const Rep& v ) const;

  // -- substraction with a scalar: substraction with I*val
  void sub   ( Rep& res, const Rep& u, const Type_t& val ) const;
  void sub   ( Rep& res, const Type_t& val, const Rep& v ) const;

  // -- map of a inplace unary operator, with operator()( Type_t& res)
  template<class OP>
  void map ( Rep& res, OP& op ) const;

  // -- map of a unary operator, with operator()( Type_t& res, const Type_t& val)
  template<class OP>
  void map ( Rep& res, OP& op, const Rep& u ) const;

  // -- with operator()( Type_t& res, const Type_t& v1, const Type_t& v2)
  template<class OP>
  void map ( Rep& res, OP& op, const Rep& u, const Rep& u ) const;

  // -- IO
  istream& read ( istream& s );
  ostream& write( ostream& s ) const;
  istream& read ( istream& s, Rep& r ) const;
  ostream& write( ostream& s, const Rep& r ) const;
};

#endif
