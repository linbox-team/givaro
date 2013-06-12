// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/matrix/givmatsparseops.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givmatsparseops.inl,v 1.3 2011-02-02 16:23:56 bboyer Exp $
// ==========================================================================

#error "this looks very much like dead code"

namespace Givaro {
#pragma message "#warning this file will probably not compile"

  // -- map of a unary operator, with operator()( Type_t& res )
  // res and u could be aliases if OP permits it
template<class Domain>
template<class UNOP>
inline void MatrixDom<Domain,Sparse>::
  map ( Rep& res, UNOP& op ) const
{
  size_t sz = res._data.size();
  for (size_t i=0; i<sz; ++i) op(res._data[i]);
}


template<class Domain>
template<class UNOP>
inline void MatrixDom<Domain,Sparse>::
  map ( Rep& res, UNOP& op, const Rep& u ) const
{
  res.reallocate(nrow(u), ncol(u));
  res._index.reallocate(u._index.size());
  res._data.reallocate(u._data.size());
  size_t sz = res._data.size();
  for (size_t i=0; i<sz; ++i) {
    res._index[i] = u._index[i];
    op(res._data[i], u._data[i]);
  }
}

  // -- Comparaizon
template<class Domain>
inline int MatrixDom<Domain,Sparse>::areEqual
 ( const Rep& P, const Rep& Q) const
{
  if (nrow(P) != nrow(Q)) return 0;
  if (ncol(P) != ncol(Q)) return 0;
  size_t sz = P._data.size();
  if (sz != Q._data.size()) return 0;
  for(size_t i=0; i<sz; ++i) {
    if (P._index[i] != Q._index[i]) return 0;
    if (_domain.areNEqual(P._data[i] != Q._data[i])) return 0;
  }
  return 1;
}

template<class Domain>
inline int MatrixDom<Domain,Sparse>::areNEqual
 ( const Rep& P, const Rep& Q) const
{
  return !areEqual(P,Q);
}

template<class Domain>
inline int MatrixDom<Domain,Sparse>::iszero ( const Rep& P ) const
{
  if (nrow(P) == 0) return 1; // -- col souhld be 0
  if (ncol(P) == 0) { cerr << " Error: inconsistent data structure";
                      return 1; } // -- row souhld be 0
  size_t sz = P._data.size();
  if (sz == 0) return 1;
  for(size_t i=0; i<sz; ++i)
    if (!_domain.iszero(P._data[i])) return 0;
  return 1;
}

template<class Domain>
inline void MatrixDom<Domain,Sparse>::mulin
 ( Rep& res, const Type_t& u ) const
{
  Curried2<MulOp<Domain> > op(_domain, u);
  map(res, op);
}

template<class Domain>
inline void MatrixDom<Domain,Sparse>::mul
 ( Rep& res, const Type_t& u, const Rep& v ) const
{
  Curried1<MulOp<Domain> > op(_domain, u);
  map(res, op, v);
}

template<class Domain>
inline void MatrixDom<Domain,Sparse>::mul
 ( Rep& res, const Rep& u, const Type_t& v ) const
{
  Curried2<MulOp<Domain> > op(_domain, v);
  map(res, op, u);
}


template<class Domain>
inline void MatrixDom<Domain,Sparse>::neg ( Rep& R, const Rep& P ) const
{
  size_t sz = P._data.size();
  if (R._data.size() != sz) {
    R.reallocate(nrow(P), ncol(P));
    R._index.copy(P._index);
    R._data.reallocate(P._data);
  }
  for(size_t i=0; i<sz; ++i) _domain.neg(R._data[i], P._data[i]);
}

  // VD is the vector domain for res and u
template<class Domain>
inline void MatrixDom<Domain,Sparse>::mul
 ( VectorDom<Domain,Dense>::Rep& R,
   const Rep& M,
   const VectorDom<Domain,Dense>& VD,
   const VectorDom<Domain,Dense>::Rep& U ) const
{
  Indice_t i,j;
  Indice_t irow, erow;
  for (i = nrow(M); --i>=0;)
  {
    // -- update the i-th row of R
    _domain.assign(R[i], _domain.zero);
    irow = M._rows[i];   // - first index of the ith row
    erow = M._rows[i+1]; // - last index of the ith row
    for (; irow != erow; ++irow)
      _domain.axpyin(R[i], M._data[irow], U[M._index[irow]]);
  }
}

template<class Domain>
inline void MatrixDom<Domain,Sparse>::multrans
 ( typename VectorDom<Domain,Dense>::Rep& R,
   const Rep& M,
   const VectorDom<Domain,Dense>& VS,
   const typename VectorDom<Domain,Dense>::Rep& U ) const
{
  Indice_t i,j;
  Indice_t irow, erow;
  for (i = 0; i<ncol(M); ++i)
    _domain.assign(R[i], _domain.zero);
  for (i = nrow(M); --i>=0;)
  {
    // -- update the i-th row of R
    irow = M._rows[i];   // - first index of the ith row
    erow = M._rows[i+1]; // - last index of the ith row
    for (; irow != erow; ++irow)
      _domain.axpyin(R[M._index[irow]], U[i], M._data[irow]);
  }
}



template<class Domain>
void MatrixDom<Domain, Sparse>::compact
 ( Rep& Ms,
   const MatrixDom<Domain, Dense>& MD,
   const MatrixDom<Domain, Dense>::Rep& Md)
{
  // -- Should compare _domain and MD.subdomain(): to be equal!

  // -- Iterate by row M and store non nul entry
  Indice_t next_rows =0; // next entry to write in _storage._rows
  Indice_t next_val = 0; // next entry to write in _data & _index
  Indice_t nrows = MD.nrow(Md);
  Indice_t ncols = MD.ncol(Md);
  size_t size = 0;   // size of _data and _index
  Ms.reallocate( nrows, ncols );
  Indice_t i,j;
  Ms._rows[0] = 0;
  Type_t val;
  for (i=0; i<nrows; ++i)
  {
    Ms._rows[1+i] = Ms._rows[i];
    for (j=0; j<ncols; ++j)
    {
      _domain.assign(val, Md(i,j));
      if ( !_domain.iszero(val) )
      {
        Ms._data.reallocate(size+1);
        Ms._index.reallocate(size+1);
        _domain.assign(Ms._data[size], val);
        Ms._index[size] = j;
        ++size; Ms._rows[i+1] +=1;
      }
    }

    cout << "i:" << i << "----> ";
    for (Indice_t k=0; k<=nrows; ++k)
      cout << "," << Ms._rows[k];
    cout << endl;
  }
}

template <class Domain>
inline ostream& MatrixDom<Domain,Sparse>::write( ostream& sout ) const
{
  return _domain.write(sout << '(') << ",Sparse)";
}

template <class Domain>
inline istream& MatrixDom<Domain,Sparse>::read( istream& sin )
{
  char ch;
  sin >> std::ws >> ch;
  if (ch != '(')
    GivError::throw_error(
      GivBadFormat("MatrixDom<Domain,Sparse>::read: syntax error no '('"));

  // -- read subdomain:
  _domain.read(sin);

  // -- read ,
  sin >> std::ws >> ch;
  if (ch != ',')
    GivError::throw_error(
      GivBadFormat("MatrixDom<Domain,Sparse>::read: syntax error no ','"));

  // -- read dense:
  char name[7];
  sin >> std::ws; sin.read(name, 6); name[6] = (char)0;
  if (strcmp(name, "Sparse") !=0)
    GivError::throw_error(
      GivBadFormat("MatrixDom<Domain,Sparse>::read: syntax error no 'Sparse'"));

  // -- read )
  sin >> std::ws >> ch;
  if (ch != ')')
    GivError::throw_error(
      GivBadFormat("MatrixDom<Domain,Sparse>::read: syntax error no ')'"));
  return sin;
}

template <class Domain>
ostream& MatrixDom<Domain,Sparse>::write( ostream& sout, const Rep& A) const
{
  sout << '[' << nrow(A) << ',' << ncol(A) << ',';
  IntDom D;
{
  VectorDom<IntDom,Dense> VDi (D);
  VDi.write(sout, A._rows)  << ',';
  VDi.write(sout, A._index) << ',';
}{
  const VectorDom<Domain,Dense> VD (_domain);
  VD.write(sout, A._data);
}
  return sout << ']';
}

template <class Domain>
istream& MatrixDom<Domain,Sparse>::read( istream& sin, Rep& R) const
{
  long nr,nc;
  char ch;
  sin >> std::ws >> ch;
  if (ch != '[')
    GivError::throw_error(
      GivBadFormat("MatrixDom<Domain,Sparse>::read: syntax error no '['"));

  sin >> std::ws >> nr;
  sin >> std::ws >> ch;
  if (ch != ',')
    GivError::throw_error(
      GivBadFormat("MatrixDom<Domain,Sparse>::read: syntax error no ','"));
  sin >> std::ws >> nc;

  R.reallocate(nr,nc);

  sin >> std::ws >> ch;
  if (ch != ',')
    GivError::throw_error(
      GivBadFormat("MatrixDom<Domain,Sparse>::read: syntax error no ','"));

  VectorDom<IntDom,Dense> VDi = IntDom();
  VectorDom<Domain,Dense> VD (_domain);
  VDi.read(sin, R._rows)  >> std::ws >> ch;
  if (ch != ',')
    GivError::throw_error(
      GivBadFormat("MatrixDom<Domain,Sparse>::read: syntax error no ','"));
  VDi.read(sin, R._index) >> std::ws >> ch;
  if (ch != ',')
    GivError::throw_error(
      GivBadFormat("MatrixDom<Domain,Sparse>::read: syntax error no ','"));
  VD.read(sin, R._data) >> std::ws >> ch;
  if (ch != ']')
    GivError::throw_error(
      GivBadFormat("MatrixDom<Domain,Sparse>::read: syntax error no ']'"));
  return sin;
}

} // Givaro

