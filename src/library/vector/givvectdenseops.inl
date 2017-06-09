// ==========================================================================
// $Source
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id
// ==========================================================================
namespace Givaro {
#pragma message "#warning this file will probably not compile"


  // -- map of a unary operator, with operator()( Type_t& res )
  // res and u could be aliases if OP permits it
template<class Domain>
template<class UNOP>
inline void VectorDom<Domain,Dense>::
  map ( Rep& res, UNOP& op ) const
{
  size_t sz = dim(res);
  for (size_t i=0; i<sz; ++i) op(res[i]);
}


template<class Domain>
template<class BINOP>
inline void VectorDom<Domain,Dense>::
  map( Rep& res, const BINOP& OP, const Rep& u, const Rep& v ) const
{
  GIVARO_ASSERT((dim(res) == dim(u)) && (dim(u) == dim(v)), "Bad size");
  size_t sz = dim(res);
  for (size_t i=0; i<sz; ++i)
    OP(res[i], u[i], v[i]);
}

template<class Domain>
template<class UNOP>
inline void VectorDom<Domain,Dense>::
  map( Rep& res, UNOP& OP, const Rep& u ) const
{
  GIVARO_ASSERT((dim(res) == dim(u)), "Bad size");
  size_t sz = dim(u);
  for (size_t i=0; i<sz; ++i)
    OP(res[i], u[i]);
}

// --- mul
template<class Domain>
inline void VectorDom<Domain,Dense>::mulin( Rep& res, const Type_t& val ) const
{
  GIVARO_ASSERT((dim(res) == dim(u)), "Bad size");
  Curried2<MulOp<Domain> > opcode(_domain, (Type_t&)val);
  map( res, opcode);
}

template<class Domain>
inline void VectorDom<Domain,Dense>::mul
 ( Rep& res, const Rep& u, const Type_t& val ) const
{
  GIVARO_ASSERT((dim(res) == dim(u)), "Bad size");
  Curried2<MulOp<Domain> > opcode(_domain, (Type_t&)val);
  map( res, opcode, u);
}

template<class Domain>
inline void VectorDom<Domain,Dense>::mul
 ( Rep& res, const Type_t& val, const Rep& v ) const
{
  GIVARO_ASSERT((dim(res) == dim(v)), "Bad size");
  Curried1<MulOp<Domain> > opcode(_domain, val);
  map( res, opcode, v);
}



// --- add
template<class Domain>
inline void VectorDom<Domain,Dense>::addin( Rep& res, const Rep& u ) const
{
  GIVARO_ASSERT((dim(res) == dim(u)), "Bad size");
  AddOp<Domain> opcode( _domain );
  map( res, opcode, u );
}

template<class Domain>
inline void VectorDom<Domain,Dense>::add
 ( Rep& res, const Rep& u, const Rep& v ) const
{
  GIVARO_ASSERT((dim(res) == dim(u)) && (dim(u) == dim(v)), "Bad size");
  AddOp<Domain> opcode (_domain);
  map( res, opcode, u, v);
}

template<class Domain>
inline void VectorDom<Domain,Dense>::add
 ( Rep& res, const Rep& u, const Type_t& val ) const
{
  GIVARO_ASSERT((dim(res) == dim(u)), "Bad size");
  Curried2<AddOp<Domain> > opcode(_domain, (Type_t&)val);
  map( res, opcode, u);
}

template<class Domain>
inline void VectorDom<Domain,Dense>::add
 ( Rep& res, const Type_t& val, const Rep& v ) const
{
  GIVARO_ASSERT((dim(res) == dim(v)), "Bad size");
  Curried1<AddOp<Domain> > opcode(_domain, val);
  map( res, opcode, v);
}

template<class Domain>
inline void VectorDom<Domain,Dense>::add
 ( Rep& res, const VectorDom<Domain,Sparse>::Rep& u, const Rep& v ) const
{
  GIVARO_ASSERT((dim(res) == dim(v)) && (dim(v) == u.dim()), "Bad size");
  // -- here assume : subdomains of res, u and v are equal
  assign(res, v);
  addin(res, u);
}

template<class Domain>
inline void VectorDom<Domain,Dense>::add
( Rep& res, const Rep& v, const VectorDom<Domain,Sparse>::Rep& u ) const
{
  GIVARO_ASSERT((dim(res) == dim(v)) && (dim(v) == u.dim()), "Bad size");
  // -- here assume : subdomains of res, u and v are equal
  assign(res, v);
  addin(res,u);
}

template<class Domain>
inline void VectorDom<Domain,Dense>::addin
 ( Rep& res, const VectorDom<Domain,Sparse>::Rep& u ) const
{
  GIVARO_ASSERT((dim(res) == u.dim()), "Bad size");
  // -- here assume : subdomains of res, u and v are equal
  typedef VectorDom<Domain,Sparse> VSparseDom;
  typedef VSparseDom::constIterator_t cIterator_t;
  typedef VSparseDom::IndiceIterator_t cIndiceIterator_t;
  cIterator_t u_curr = u.begin_data();
  cIterator_t u_end = u.end_data();
  cIndiceIterator_t u_indice = u.begin_indice();
  for ( ; u_curr != u_end; ++u_curr, ++u_indice )
    _domain.addin(res[*u_indice], *u_curr);
}



// --- sub
template<class Domain>
inline void VectorDom<Domain,Dense>::subin( Rep& res, const Rep& u ) const
{
  GIVARO_ASSERT((dim(res) == dim(u)), "Bad size");
  SubOp<Domain> opcode ( _domain );
  map( res, opcode, u );
}

template<class Domain>
inline void VectorDom<Domain,Dense>::sub( Rep& res, const Rep& u, const Rep& v ) const
{
  GIVARO_ASSERT((dim(res) == dim(v)) && (dim(v) == dim(u)), "Bad size");
  SubOp<Domain> opcode ( _domain );
  map( res, opcode, u, v);
}

template<class Domain>
inline void VectorDom<Domain,Dense>::sub( Rep& res, const Rep& u, const Type_t& val ) const
{
  GIVARO_ASSERT((dim(res) == dim(u)), "Bad size");
  Curried2<SubOp<Domain> > opcode( _domain, (Type_t&)val);
  map( res, opcode, u);
}

template<class Domain>
inline void VectorDom<Domain,Dense>::sub( Rep& res, const Type_t& val, const Rep& v ) const
{
  GIVARO_ASSERT((dim(res) == dim(v)), "Bad size");
  Curried1<SubOp<Domain> > opcode( _domain, val);
  map( res, opcode, v);
}

template<class Domain>
inline void VectorDom<Domain,Dense>::subin
 ( Rep& res, const VectorDom<Domain,Sparse>::Rep& u ) const
{
  GIVARO_ASSERT((dim(res) == u.dim()), "Bad size");
  // -- here assume : subdomains of res, u and v are equal
  typedef VectorDom<Domain,Sparse> VSparseDom;
  typedef VSparseDom::constIterator_t cIterator_t;
  typedef VSparseDom::IndiceIterator_t cIndiceIterator_t;
  cIterator_t u_curr = u.begin_data();
  cIterator_t u_end = u.end_data();
  cIndiceIterator_t u_indice = u.begin_indice();
  for ( ; u_curr != u_end; ++u_curr, ++u_indice )
    _domain.subin(res[*u_indice], *u_curr);
}

template<class Domain>
inline void VectorDom<Domain,Dense>::sub
 ( Rep& res, const VectorDom<Domain,Sparse>::Rep& u, const Rep& v ) const
{
  GIVARO_ASSERT((dim(res) == dim(v)) && (dim(v) == u.dim()), "Bad size");
  // -- here assume : subdomains of res, u and v are equal
  neg(res, v);
  addin(res, u);
}

template<class Domain>
inline void VectorDom<Domain,Dense>::sub
 ( Rep& res, const Rep& v, const VectorDom<Domain,Sparse>::Rep& u ) const
{
  GIVARO_ASSERT((dim(res) == dim(v)) && (dim(v) == u.dim()), "Bad size");
  // -- here assume : subdomains of res, u and v are equal
  assign(res,v);
  subin(res,u);
}


// --- sub

// --- neg
template<class Domain>
inline void VectorDom<Domain,Dense>::negin( Rep& res ) const
{
  NegOp<Domain> opcode ( _domain );
  map( res, opcode, res );
}

template<class Domain>
inline void VectorDom<Domain,Dense>::neg( Rep& res, const Rep& u ) const
{
  GIVARO_ASSERT((dim(res) == dim(u)), "Bad size");
  NegOp<Domain> opcode ( _domain );
  map( res, opcode, u );
}

template<class Domain>
void VectorDom<Domain,Dense>::axpy
 ( Rep& res, const Type_t& a, const Rep& x, const Rep& y ) const
{
  GIVARO_ASSERT((dim(res) == dim(x)) && (dim(x) == dim(y)), "Bad size");
  size_t sz = dim(res);
  for ( size_t i=0; i<sz; ++i )
   _domain.axpy(res[i], a, x[i], y[i]);
}

template<class Domain>
void VectorDom<Domain,Dense>::axpyin
 ( Rep& res, const Type_t& a, const Rep& x ) const
{
  GIVARO_ASSERT((dim(res) == dim(x)), "Bad size");
  size_t sz = dim(res);
  for ( size_t i=0; i<sz; ++i )
   _domain.axpyin(res[i], a, x[i]);
}

template<class Domain>
void VectorDom<Domain,Dense>::axmy
 ( Rep& res, const Type_t& a, const Rep& x, const Rep& y ) const
{
  GIVARO_ASSERT((dim(res) == dim(x)) && (dim(x) == dim(y)), "Bad size");
  size_t sz = dim(res);
  for ( size_t i=0; i<sz; ++i )
   _domain.axmy(res[i], a, x[i], y[i]);
}

template<class Domain>
void VectorDom<Domain,Dense>::axmyin
 ( Rep& res, const Type_t& a, const Rep& x ) const
{
  GIVARO_ASSERT((dim(res) == dim(x)), "Bad size");
  size_t sz = dim(res);
  for ( size_t i=0; i<sz; ++i )
   _domain.axmyin(res[i], a, x[i]);
}


template<class Domain>
void VectorDom<Domain,Dense>::axpy
 ( Rep& res, const Rep& a, const Rep& b, const Rep& y ) const
{
  GIVARO_ASSERT((dim(res) == dim(a)) && (dim(a) == dim(b)) && (dim(b) == dim(y)), "Bad size");
  size_t sz = dim(res);
  for ( size_t i=0; i<sz; ++i )
   _domain.axpy(res[i], a[i], b[i], y[i]);
}

template<class Domain>
void VectorDom<Domain,Dense>::axpyin
 ( Rep& res, const Rep& a, const Rep& b ) const
{
  GIVARO_ASSERT((dim(res) == dim(a)) && (dim(a) == dim(b)), "Bad size");
  size_t sz = dim(res);
  for ( size_t i=0; i<sz; ++i )
   _domain.axpyin(res[i], a[i], b[i]);
}

template<class Domain>
void VectorDom<Domain,Dense>::axmy
 ( Rep& res, const Rep& a, const Rep& b, const Rep& y ) const
{
  GIVARO_ASSERT((dim(res) == dim(a)) && (dim(a) == dim(b)) && (dim(b) == dim(y)), "Bad size");
  size_t sz = dim(res);
  for ( size_t i=0; i<sz; ++i )
   _domain.axmy(res[i], a[i], b[i], y[i]);
}

template<class Domain>
void VectorDom<Domain,Dense>::axmyin
 ( Rep& res, const Rep& a, const Rep& b ) const
{
  GIVARO_ASSERT((dim(res) == dim(a)) && (dim(a) == dim(b)), "Bad size");
  size_t sz = dim(res);
  for ( size_t i=0; i<sz; ++i )
   _domain.axmyin(res[i], a[i], b[i]);
}



template<class Domain>
inline void VectorDom<Domain,Dense>::dot
  ( Type_t& res, const Rep& op1, const Rep& op2) const
{
  GIVARO_ASSERT((dim(op1) == dim(op2)), "Bad size");
  size_t sz = dim(op1);
  _domain.assign(res, _domain.zero);
  for ( size_t i=0; i<sz; ++i ) _domain.axpyin(res, op1[i], op2[i] );
}



// ==========================================================================
//
// -- Write the domain
//
template<class Domain>
ostream& VectorDom<Domain, Dense>::write( ostream& o ) const
{
  return _domain.write(o << "(") << ",Dense)";
}

// -- read the domain
template<class Domain>
istream& VectorDom<Domain, Dense>::read( istream& sin )
{
  char ch;
  sin >> std::ws >> ch;
  if (ch != '(')
    GivError::throw_error(
      GivBadFormat("VectorDom<Domain,Dense>::read: syntax error no '('"));

  // -- read subdomain:
  _domain.read(sin);

  // -- read ,
  sin >> std::ws >> ch;
  if (ch != ',')
    GivError::throw_error(
      GivBadFormat("VectorDom<Domain,Dense>::read: syntax error no ','"));

  // -- read dense:
  char name[6];
  sin >> std::ws; sin.read(name, 5); name[5] = (char)0;
  if (strcmp(name, "Dense") !=0)
    GivError::throw_error(
      GivBadFormat("VectorDom<Domain,Dense>::read: syntax error no 'Dense'"));

  // -- read )
  sin >> std::ws >> ch;
  if (ch != ')')
    GivError::throw_error(
      GivBadFormat("VectorDom<Domain,Dense>::read: syntax error no ')'"));
  return sin;
}


// ==========================================================================
//
// -- Don't write the domain !, only write a rep / domain
template<class Domain>
ostream& VectorDom<Domain, Dense>::write(ostream& o, const Rep& V) const
{
   if (dim(V) ==0) return o << "[ ]";
   o << '[';
   size_t i;
   for (i=0; i<dim(V)-1; i++) o << V[i] << ',';
   return o << V[i] << ']';
}


//
// Read a vector given by a Maple-like format:
// the grammar is :
//   s  --->  '[' list_of_elt ']'
//   list_of_elt --->   elt
//                    | list_of_elt ',' elt
// The contraints are :
//   All lines are the same number of Elements.
//   The separators a those of the C lexical-convention, i.e.
//    ' ', '\n', '\t', '\f' .

template<class Domain>
istream&  VectorDom<Domain,Dense>::read (istream& fin, Rep& A) const
{
  char ch;

  // -- Skip the first "white":
  fin >> std::ws; fin.get(ch);
  if (ch != '[')
    GivError::throw_error(
      GivBadFormat("VectorDom<Domain,Dense>::read: syntax error no '['"));

  // -- Read the line and count the nb of elts
  int i = 0;
  Rep rep;
  rep.allocate(1);
  Type_t Tmp;
  fin >> Tmp;
  fin >> std::ws; fin.get(ch);
  rep[0] = Tmp;
  while (ch != ']')
  {
    if (ch != ',')
      GivError::throw_error(
        GivBadFormat("VectorDom<T,Dense>::read: syntax error no ','"));
    i++;
    fin >> std::ws >> Tmp;
    fin >> std::ws; fin.get(ch);

    // resize the vector :
    rep.resize(i+1);
    rep[i] = Tmp;
  }
  A.logcopy( rep );
  return fin;
}
} // Givaro
#include "givaro/givvectdensespe.inl"
