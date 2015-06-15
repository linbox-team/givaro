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
#pragma message "#warning this file will probably not compile"



template<class Domain>
class VectorDom<Domain, Sparse> {
  Domain _domain;	// domain of the entries
public :
  // -- Exported types
  typedef typename Domain::Rep   			 	Type_t;
  typedef 	   Domain   				 	Domain_t;
  typedef 	   int 					 	Indice_t;
  typedef 	   Sparse	 			 	StorageTag_t;
  typedef typename RetVectorStorage<Type_t,Sparse>::Storage_t 	Storage_t;

  // -- Representation of Element of VectorDom<D, Sparse>
  typedef 	   Storage_t 	 				Rep;

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
  void init ( Rep& r, size_t dim =0) const
  { r.allocate(dim,0); }

  // -- assignment operator: from a vector
  void assign (Rep& r, const Rep& v)
  {
    r.copy(v);
  }

  // -- Comparaizon
  int areEqual ( const Rep& P, const Rep& Q) const;
  int areNEqual( const Rep& P, const Rep& Q) const;
  int iszero  ( const Rep& P ) const;

  // -- return the dimension of a vector
  size_t dim( const Rep& u ) const { return u.size(); }
  const Domain& subdomain() const { return _domain; }

  // -- Arithmetic operations: base
  void add ( Rep& res, const Rep& op1, const Rep& op2) const;
  void sub ( Rep& res, const Rep& op1, const Rep& op2) const;

  // -- dot product: operands could be aliased
  void dot ( Type_t& res, const Rep& u, const Rep& v ) const;

  // -- Syntaxic sugar: (Value) op (Vector): Element wise ops.
  void addin( Rep& res, const Rep& u ) const;
  void add  ( Rep& res, const Rep& u, const Type_t& val ) const;
  void add  ( Rep& res, const Type_t& val, const Rep& v ) const;
  void subin( Rep& res, const Rep& u ) const;
  void sub  ( Rep& res, const Rep& u, const Type_t& val ) const;
  void sub  ( Rep& res, const Type_t& val, const Rep& v ) const;
  void negin( Rep& res ) const;
  void neg  ( Rep& res, const Rep& u ) const;

  // -- Compression method to compact a dense vector
  void compact( Rep& u, const VectorDom<Domain, Dense>& VDom,
                const typename VectorDom<Domain, Dense>::Rep& v ) const;

  // -- Compression method to compact a sparse vector
  void compact( Rep& u, const VectorDom<Domain, Sparse>& VDom,
                const typename VectorDom<Domain, Sparse>::Rep& v ) const;

  template<class UNOP>
  void map( Rep& r, const UNOP& op, const Rep& u) const;

  template<class UNOP>
  void map( Rep& r, UNOP& op, const Rep& u) const;

  // -- IO: domain
  ostream& write( ostream& o ) const;
  istream& read ( istream& i );

  // -- IO: domain Element
  ostream& write( ostream& o, const Rep& r ) const;
  istream& read ( istream& i, Rep& r ) const;


  // -- Iteration over a sparse vector:
  typedef typename RetVectorStorage<Type_t,Sparse>::Iterator_t 		Iterator_t;
  typedef typename RetVectorStorage<Type_t,Sparse>::constIterator_t 	constIterator_t;
  typedef typename RetVectorStorage<Type_t,Sparse>::IndiceIterator_t 	IndiceIterator_t;

  Iterator_t	  begin_data( Rep& U ) const { return U.begin_data(); }
  Iterator_t	  end_data  ( Rep& U ) const { return U.end_data(); }
  constIterator_t begin_data( const Rep& U ) const { return U.begin_data(); }
  constIterator_t end_data  ( const Rep& U ) const { return U.end_data(); }
  IndiceIterator_t begin_indice( const Rep& U ) const { return U.begin_indice(); }
  IndiceIterator_t end_indice  ( const Rep& U ) const { return U.end_indice(); }
};

} //Givaro
#endif
