#ifndef _GIV_OPERATION_H_
#define _GIV_OPERATION_H_
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/tools/givops.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givops.h,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description:
// Class version of operations on a group F. By default use arithmetic operator
// 
#include "givaro/givconfig.h"
#include "givaro/giviterator.h"

template<class Domain>
struct BaseOP {
  typedef typename 	Domain::Rep 	Type_t;
  typedef 		Domain 		Domain_t;
  const Domain& _domain;
  BaseOP( const Domain& D ) : _domain(D) {} 
};

template<class Domain>
struct CopyOp : public BaseOP<Domain> {
  typedef typename 	BaseOP<Domain>::Type_t 	Type_t;
  typedef 		Domain 			Domain_t;
  CopyOp( const Domain& D ) : BaseOP<Domain>(D) {}
  void operator() (Type_t& v1, const Type_t& v2) const {  
    _domain.assign(v1, v2); 
  }
};

// -- Default operators

template<class Domain>
struct MulOp : public BaseOP<Domain> {
  typedef typename 	BaseOP<Domain>::Type_t 	Type_t;
  typedef 		Domain 			Domain_t;
  MulOp( const Domain& D ) : BaseOP<Domain>(D) {}
  void operator() (Type_t& res, const Type_t& v1, const Type_t& v2) const 
  { _domain.mul( res, v1, v2); }
  void operator() (Type_t& res, const Type_t& v1 ) const 
  { _domain.mulin( res, v1 ); }
};

template<class Domain>
struct DivOp : public BaseOP<Domain> {
  typedef typename 	BaseOP<Domain>::Type_t 	Type_t;
  typedef 		Domain 			Domain_t;
  DivOp( const Domain& D ) : BaseOP<Domain>(D) {}
  void operator() (Type_t& res, const Type_t& v1, const Type_t& v2) const 
  { _domain.div( res, v1, v2); }
  void operator() (Type_t& res, const Type_t& v1) const 
  { _domain.divin( res, v1 ); }
};

template<class Domain>
struct ModOp : public BaseOP<Domain> {
  typedef typename 	BaseOP<Domain>::Type_t 	Type_t;
  typedef 		Domain 			Domain_t;
  ModOp( const Domain& D ) : BaseOP<Domain>(D) {}
  void operator() (Type_t& res, const Type_t& v1, const Type_t& v2) const 
  { _domain.mod( res, v1, v2); }
  void operator() (Type_t& res, const Type_t& v1) const 
  { _domain.modin( res, v1 ); }
};

template<class Domain>
struct AddOp : public BaseOP<Domain> {
  typedef typename 	BaseOP<Domain>::Type_t 	Type_t;
  typedef 		Domain 			Domain_t;
  AddOp( const Domain& D ) : BaseOP<Domain>(D) {}
  void operator() (Type_t& res, const Type_t& v1, const Type_t& v2) const 
  { _domain.add( res, v1, v2 ); }
  void operator() (Type_t& res, const Type_t& v1) const 
  { _domain.addin( res, v1 ); }
};

template<class Domain>
struct SubOp : public BaseOP<Domain> {
  typedef typename 	BaseOP<Domain>::Type_t 	Type_t;
  typedef 		Domain 			Domain_t;
  SubOp( const Domain& D ) : BaseOP<Domain>(D) {}
  void operator() (Type_t& res, const Type_t& v1, const Type_t& v2) const 
  { _domain.sub( res, v1, v2 ); }
  void operator() (Type_t& res, const Type_t& v1) const 
  { _domain.subin( res, v1 ); }
};

template<class Domain>
struct NegOp : public BaseOP<Domain> {
  typedef typename 	BaseOP<Domain>::Type_t 	Type_t;
  typedef 		Domain 			Domain_t;
  NegOp( const Domain& D ) : BaseOP<Domain>(D) {}
  void operator() (Type_t& res, const Type_t& v1) const 
  { _domain.neg( res, v1 ); }
  void operator() (Type_t& res ) const 
  { _domain.negin( res, res ); }
};

template<class Domain>
struct MulAddOp : public BaseOP<Domain> {
  typedef typename 	BaseOP<Domain>::Type_t 	Type_t;
  typedef 		Domain 			Domain_t;
  MulAddOp( const Domain& D ) : BaseOP<Domain>(D) {}
  void operator()(Type_t& res, const Type_t& v1, const Type_t& v2, const Type_t& v3) const
  { _domain.axpy( res, v1, v2, v3 ); }
  void operator() (Type_t& res, const Type_t& v1, const Type_t& v2 ) const 
  { _domain.axpy( res, v1, v2, res ); }
};


// -- BinOP -> UnOP
// -- UnOp -> Zero OP
template<class OP>
struct Curried1 : public OP {
  typedef typename OP::Type_t   Type_t;
  typedef typename OP::Domain_t Domain_t;
  Type_t& _val;
  Curried1( const Domain_t& D, Type_t& val ) : OP(D), _val(val) {}
  Curried1( const Domain_t& D, const Type_t& val ) : OP(D), _val((Type_t&)val) {}
  void operator()(Type_t& v1, const Type_t& v2) { OP::operator()(v1, _val, v2); }
  void operator()(Type_t& v1) { OP::operator()(v1, _val); }
};

// -- BinOP -> UnOP
template<class OP>
struct Curried2 : public OP {
  typedef typename OP::Type_t   Type_t;
  typedef typename OP::Domain_t Domain_t;
  const Type_t& _val;
  Curried2( const Domain_t& D, Type_t& val ) : OP(D), _val(val) {}
  Curried2( const Domain_t& D, const Type_t& val ) : OP(D), _val((Type_t&)val) {}
  void operator()(Type_t& v1, const Type_t& v2) { OP::operator()(v1, v2, _val); }
  void operator()(Type_t& v1) { OP::operator()(v1, _val); }
};


#endif
