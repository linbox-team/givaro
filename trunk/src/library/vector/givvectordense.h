#ifndef _GIV_VECTOR_DENSE_H_
#define _GIV_VECTOR_DENSE_H_
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/vector/givvectordense.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givvectordense.h,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description:
// Domain of dense vector over K with classic arithmetic operations
// over T (vector x vector, vector x T, scalar product, shift).
// A element of this domain is a vector of K^n for any n.
// Vector handle computation over sub part of continuous elements of
// a vector as well as stride.

#include "givaro/givvector.h"
#include "givaro/givvectorsparse.h"
#include "givaro/givstoragedense.h"
#include "givaro/givelem.h"


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


  // Vector (+/-/*) Value ==  element wise operation 
  void mulin( Rep& res, const Type_t& u ) const;
  void mul  ( Rep& res, const Rep& u, const Type_t& val ) const;
  void mul  ( Rep& res, const Type_t& val, const Rep& v ) const;

  void add ( Rep& res, const Rep& u, const Type_t& val ) const;
  void add ( Rep& res, const Type_t& val, const Rep& v ) const;

  void sub  ( Rep& res, const Rep& u, const Type_t& val ) const;
  void sub  ( Rep& res, const Type_t& val, const Rep& v ) const;


  // - axpy like operations, element wise:
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
  template<class OP>
  void map ( Rep& res, OP& op ) const;

  // -- map of a unary operator, with operator()( Type_t& res, const Type_t& val)
  // res and u could be aliases if OP permits it
  template<class OP>
  void map ( Rep& res, OP& op, const Rep& u ) const;

  // -- map of a binary operator, with :
  // -- operator()( Type_t& res, const Type_t&, const Type_t& )
  template<class OP>
  void map ( Rep& res, const OP& op, const Rep& u, const Rep& u ) const;

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



#endif
