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
#pragma message "#warning this file will probably not compile"



template<class Domain>
class VectorDom<Domain, Dense> {
public:
  Domain 	_domain; // domain of the entries
public:
  // -- Exported types
  typedef 	   Domain 					Domain_t;
  typedef typename Domain::Rep 					Type_t;
  typedef 	   int 						Indice_t;
  typedef 	   Dense	 				StorageTag_t;
  typedef typename RetVectorStorage<Type_t,Dense>::Storage_t 	Storage_t;

  // -- The representation of a dense vector
  typedef 	   Storage_t					Rep;

  // -- Self_t
  typedef 	   VectorDom<Domain, Dense> 			Self_t;

  // -- Dstor
  ~VectorDom() {}

  // -- Cstor of a new vector space of dimension s
  VectorDom<Domain,Dense>() : _domain() {}

  VectorDom<Domain,Dense>( const Domain& dom ) : _domain(dom) {}

  // -- Cstor of recopy
  VectorDom<Domain,Dense>( const Self_t& V ) : _domain(V._domain) {}

  int operator==( const VectorDom<Domain,Dense>& BC) const
  { return _domain == BC._domain;}
  int operator!=( const VectorDom<Domain,Dense>& BC) const
  { return _domain != BC._domain;}


  // -- init :
  void init( Rep& v, size_t dim =0 ) const
  {
    v.reallocate(dim);
  }
  // --
  void init( Rep& v, const Rep& u ) const
  {
    v.copy(v);
  }

  // -- assignment operator: from a vector of the same vect space
  void assign ( Rep& r, const Rep& v) const
  {
    r.copy(v);
  }
  // -- assignment operator: from value * [1,.....1]
  void assign ( Rep& r, size_t dim, const Type_t& val) const
  {
    r.reallocate( dim );
    for (size_t i=0; i<dim; ++i) r[i] = val;
  }

  // -- Comparizon
  int areEqual (const Rep& P, const Rep& Q) const
  {
    size_t sP = P.size(), sQ = Q.size();
    if (sP != sQ) return 0;
    for (int i=0; i<sP; ++i)
      if (!_domain.areEqual(P[i], Q[i])) return 0;
    return 1;
  }
  int areNEqual(const Rep& P, const Rep& Q) const
  {
    return !areEqual(P,Q);
  }
  int iszero(const Rep& P) const
  {
    size_t sP =P.size();
    if (sP ==0) return 1;
    for (int i=0; i<sP; ++i)
      if (!_domain.iszero(P[i])) return 0;
    return 1;
  }

  // -- Return the domain of the entries
  const Domain& subdomain() const { return _domain; }

  // -- return the dimension of the vector space
  size_t dim( const Rep& v ) const { return v.size(); }

  // -- arithmetic operator: operands could be aliased
  void addin ( Rep& res, const Rep& u ) const;
  void add ( Rep& res, const Rep& u, const Rep& v ) const;

  void addin( Rep& res, const VectorDom<Domain,Sparse>::Rep& v ) const;
  void add ( Rep& res, const Rep& u, const VectorDom<Domain,Sparse>::Rep& v ) const;
  void add ( Rep& res, const VectorDom<Domain,Sparse>::Rep& u, const Rep& v ) const;

  void subin( Rep& res, const Rep& u ) const;
  void sub  ( Rep& res, const Rep& u, const Rep& v ) const;

  void subin( Rep& res, const VectorDom<Domain,Sparse>::Rep& v ) const;
  void sub ( Rep& res, const Rep& u, const VectorDom<Domain,Sparse>::Rep& v ) const;
  void sub ( Rep& res, const VectorDom<Domain,Sparse>::Rep& u, const Rep& v ) const;

  void negin ( Rep& res ) const;
  void neg ( Rep& res, const Rep& u ) const;

  // - axpy like operations:
  // r <- a*x+y
  void axpy  ( Rep& res, const Type_t& a, const Rep& x, const Rep& y )const;
  // r <- r+a*x
  void axpyin( Rep& res, const Type_t& a, const Rep& x ) const;
  // r <- y-a*x
  void axmy  ( Rep& res, const Type_t& a, const Rep& x, const Rep& y ) const;
  // r <- r-a*x
  void axmyin( Rep& res, const Type_t& a, const Rep& x ) const;


  // Vector (+/-/*) Value ==  Element wise operation
  void mulin( Rep& res, const Type_t& u ) const;
  void mul  ( Rep& res, const Rep& u, const Type_t& val ) const;
  void mul  ( Rep& res, const Type_t& val, const Rep& v ) const;

  void add ( Rep& res, const Rep& u, const Type_t& val ) const;
  void add ( Rep& res, const Type_t& val, const Rep& v ) const;

  void sub  ( Rep& res, const Rep& u, const Type_t& val ) const;
  void sub  ( Rep& res, const Type_t& val, const Rep& v ) const;


  // - axpy like operations, Element wise:
  // r <- a*x+y
  void axpy  ( Rep& res, const Rep& a, const Rep& x, const Rep& y ) const;
  // r <- r+a*x
  void axpyin( Rep& res, const Rep& a, const Rep& x ) const;
  // r <- y-a*x
  void axmy  ( Rep& res, const Rep& a, const Rep& x, const Rep& y ) const;
  // r <- r-a*x
  void axmyin( Rep& res, const Rep& a, const Rep& x ) const;


  // -- dot product: operands could be aliased
  void dot ( Type_t& res, const Rep& u, const Rep& v ) const;

  // -- map of a unary operator, with operator()( Type_t& res )
  // res and u could be aliases if OP permits it
  template<class UNOP>
  void map ( Rep& res, UNOP& op ) const;

  // -- map of a unary operator, with operator()( Type_t& res, const Type_t& val)
  // res and u could be aliases if OP permits it
  template<class UNOP>
  void map ( Rep& res, UNOP& op, const Rep& u ) const;

  // -- map of a binary operator, with :
  // -- operator()( Type_t& res, const Type_t&, const Type_t& )
  template<class BINOP>
  void map ( Rep& res, const BINOP& op, const Rep& u, const Rep& u ) const;

  // -- IO
  istream& read ( istream& s );
  ostream& write( ostream& s ) const;
  istream& read ( istream& s, Rep& r ) const;
  ostream& write( ostream& s, const Rep& r ) const;

  // -- Iteration over a sparse vector:
  typedef typename
    RetVectorStorage<Type_t,Sparse>::Iterator_t          Iterator_t;
  typedef typename
    RetVectorStorage<Type_t,Sparse>::constIterator_t     constIterator_t;
  typedef typename
    RetVectorStorage<Type_t,Sparse>::IndiceIterator_t    IndiceIterator_t;

  // -- Basic iterator:
  Iterator_t begin() { return _storage.begin(); }
  Iterator_t end()   { return _storage.end(); }
  constIterator_t begin() const { return _storage.begin(); }
  constIterator_t end() const   { return _storage.end(); }

  IndiceIterator_t begin_indice() const { return IndiceIterator_t(0); }
  IndiceIterator_t end_indice() const   { return IndiceIterator_t(_dim); }
};

} // Givaro

#endif
