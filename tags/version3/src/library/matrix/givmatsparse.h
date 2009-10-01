#ifndef _GIV_MATRIX_SPARSE_H_
#define _GIV_MATRIX_SPARSE_H_
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/matrix/givmatsparse.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givmatsparse.h,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description:
// of matrix by blocks.

#include "givaro/givmatrix.h"
#include "givaro/givvector.h"
#include "givaro/givmatstoragesparse.h"

// --
// -- Matrix class: dense matrix
// --
template <class Domain>
class MatrixDom<Domain, Sparse> 
{
  Domain _domain;

public:
  typedef 	   Domain 					Domain_t;
  typedef typename Domain::Rep					Type_t;
  typedef 	   int 						Indice_t;
  typedef 	   Dense 					StorageTag_t;
  typedef typename RetMatrixStorage<Type_t,Sparse>::Storage_t 	Storage_t;

  // -- Representation of element of the domain
  typedef 	   Storage_t    				Rep;

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
  void init(Rep& r, Indice_t nr, Indice_t nc) const
  { r.allocate(nr, nc); }
  void init(Rep& r) const
  { r.allocate(0,0); }
  void init(Rep& A, const Rep& B) const
  { A.copy(B); }

  //-- Assignment operator: physical copy
  void assign (Rep& r, const Rep& a) const  
  { r.copy(a); }

  // -- Comparaizon
  int areEqual ( const Rep& P, const Rep& Q) const;
  int areNEqual( const Rep& P, const Rep& Q) const;
  int iszero  ( const Rep& P ) const;

  //-- Dimension of the matrix space
  Indice_t nrow(const Rep& A) const { return A._nrow; }
  Indice_t ncol(const Rep& A) const { return A._ncol; }
  Domain_t subdomain() const { return _domain; }

  // -- arithmetic operator: operands could be aliased
  void mulin ( Rep& res, const Type_t& u ) const;
  void mul   ( Rep& res, const Type_t& u, const Rep& v ) const;
  void mul   ( Rep& res, const Rep& u, const Type_t& v ) const;

  // VD is the vector domain for res and u
  void mul      ( VectorDom<Domain,Dense>::Rep& res, 
                  const Rep& M,
                  const VectorDom<Domain,Dense>& VD,
                  const VectorDom<Domain,Dense>::Rep& u ) const;
  void multrans ( typename VectorDom<Domain,Dense>::Rep& res, 
                  const Rep& M,
                  const VectorDom<Domain,Dense>& VS,
                  const typename VectorDom<Domain,Dense>::Rep& u ) const;


  void negin ( Rep& P ) const
  {  
    size_t sz = P._data.size();
    for(size_t i=0; i<sz; ++i) _domain.negin(P._data[i]);
  }

  void neg   ( Rep& res, const Rep& u ) const;

  // -- map of a unary operator, with operator()( Type_t& res)
  template<class OP>
  void map ( Rep& res, OP& op ) const;

  // -- map of a unary operator, with operator()( Type_t& res, const Type_t& val)
  template<class OP>
  void map ( Rep& res, OP& op, const Rep& u ) const;

  // -- IO
  istream& read ( istream& s );
  ostream& write( ostream& s ) const;
  istream& read ( istream& s, Rep& r ) const;
  ostream& write( ostream& s, const Rep& r ) const;

  // -- Compression method to compact a dense matrix to a sparse
  // template<class StorageTag>, 
  void compact( Rep& Ms,
                const MatrixDom<Domain, Dense>& MD, 
                const MatrixDom<Domain, Dense>::Rep& Md);
};

//#include "givaro/givmatsparseops.inl"

#endif