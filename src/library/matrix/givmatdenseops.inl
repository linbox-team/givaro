// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/matrix/givmatdenseops.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givmatdenseops.inl,v 1.3 2011-01-19 18:29:09 briceboyer Exp $
// ==========================================================================

#error "this looks very much like dead code"

namespace Givaro {
#pragma message "#warning this file will probably not compile"

template<class Domain>
int MatrixDom<Domain,Dense>::areEqual ( const Rep& A, const Rep& B ) const
{
	if (ncol(A) != ncol(B)) return 0;
	if (nrow(A) != nrow(B)) return 0;
	return _supportdomain.areEqual(A,B);
}

template<class Domain>
int MatrixDom<Domain,Dense>::areNEqual ( const Rep& A, const Rep& B ) const
{
	return !areEqual(A,B);
}

template<class Domain>
int MatrixDom<Domain,Dense>::iszero ( const Rep& A ) const
{
	return _supportdoamin.iszero(A);
}

// res <- A + B; aliases of operands are allowed
template<class Domain>
void MatrixDom<Domain,Dense>::add( Rep& R, const Rep& A, const Rep& B) const
{
	GIVARO_ASSERT((ncol(R) == ncol(A)) && (ncol(A) == ncol(B)), "Bad size");
	GIVARO_ASSERT((nrow(R) == nrow(A)) && (nrow(A) == nrow(B)), "Bad size");
	_supportdomain.add(R,A,B);
}

// res <- res + A; aliases of operands are allowed
template<class Domain>
void MatrixDom<Domain,Dense>::addin( Rep& R, const Rep& A ) const
{
	GIVARO_ASSERT((ncol(R) == ncol(A)), "Bad size");
	GIVARO_ASSERT((nrow(R) == nrow(A)), "Bad size");
	_supportdomain.addin(R,A);
}

// res <- A + val; aliases of operands are allowed
template<class Domain>
void MatrixDom<Domain,Dense>::add( Rep& R, const Rep& A, const Type_t& val ) const
{
	GIVARO_ASSERT((ncol(R) == ncol(A)), "Bad size");
	GIVARO_ASSERT((nrow(R) == nrow(A)), "Bad size");
	long i;
	size_t min = (ncol(R) < nrow(R) ? ncol(R) : nrow(R));
	for (i=min; --i>=0; ) _domain.add(R(i,i), A(i,i), val);
}

// res <- val + A; aliases of operands are allowed
template<class Domain>
void MatrixDom<Domain,Dense>::add( Rep& R, const Type_t& val, const Rep& A) const
{
	GIVARO_ASSERT((ncol(R) == ncol(A)), "Bad size");
	GIVARO_ASSERT((nrow(R) == nrow(A)), "Bad size");
	long i;
	size_t min = (ncol(R) < nrow(R) ? ncol(R) : nrow(R));
	for (i=min; --i>=0; ) _domain.add(R(i,i), A(i,i), val);
}



// res <- A - B; aliases of operands are allowed
template<class Domain>
void MatrixDom<Domain,Dense>::sub( Rep& R, const Rep& A, const Rep& B) const
{
	GIVARO_ASSERT((ncol(R) == ncol(A)) && (ncol(A) == ncol(B)), "Bad size");
	GIVARO_ASSERT((nrow(R) == nrow(A)) && (nrow(A) == nrow(B)), "Bad size");
	_supportdomain.sub(R,A,B);
}

// res <- res - A; aliases of operands are allowed
template<class Domain>
void MatrixDom<Domain,Dense>::subin( Rep& R, const Rep& A ) const
{
	GIVARO_ASSERT((ncol(R) == ncol(A)), "Bad size");
	GIVARO_ASSERT((nrow(R) == nrow(A)), "Bad size");
	_supportdomain.subin(R,A);
}

// res <- res - A; aliases of operands are allowed
template<class Domain>
void MatrixDom<Domain,Dense>::sub( Rep& R, const Rep& A, const Type_t& val ) const
{
	GIVARO_ASSERT((ncol(R) == ncol(A)), "Bad size");
	GIVARO_ASSERT((nrow(R) == nrow(A)), "Bad size");
	long i;
	size_t min = (ncol(R) < nrow(R) ? ncol(R) : nrow(R));
	for (i=min; --i>=0; ) _domain.sub(R(i,i), A(i,i), val);
}

// res <- res - A; aliases of operands are allowed
template<class Domain>
void MatrixDom<Domain,Dense>::sub( Rep& R, const Type_t& val, const Rep& A) const
{
	GIVARO_ASSERT((ncol(R) == ncol(A)), "Bad size");
	GIVARO_ASSERT((nrow(R) == nrow(A)), "Bad size");
	long i;
	size_t min = (ncol(R) < nrow(R) ? ncol(R) : nrow(R));
	for (i=min; --i>=0; ) _domain.sub(R(i,i), val, A(i,i));
};

// res <- - A; aliases of operands are allowed
template<class Domain>
void MatrixDom<Domain,Dense>::neg( Rep& R, const Rep& A) const
{
	GIVARO_ASSERT((ncol(R) == ncol(A)), "Bad size");
	GIVARO_ASSERT((nrow(R) == nrow(A)), "Bad size");
	_supportdomain.neg(R,A);
}

template<class Domain>
void MatrixDom<Domain,Dense>::negin( Rep& R ) const
{
	_supportdomain.negin(R);
}


// R <- A * B, R cannot be an alias for A or B
template<class Domain>
void MatrixDom<Domain,Dense>::mul( Rep& A, const Rep& B, const Rep& C) const
{
	GIVARO_ASSERT((ncol(A) == ncol(B)) && (ncol(B) == ncol(C)), "Bad size");
	GIVARO_ASSERT((nrow(A) == nrow(B)) && (nrow(B) == nrow(C)), "Bad size");
	Indice_t k, i, j;
	int startBi, startCj, startAi;
	Indice_t nrA = nrow(A);
	Indice_t ncA = ncol(A);
	Indice_t ncB = ncol(B);

	// -- the algorithm used the fact that matrices are stored by row
	Type_t tmp;
	for (i=0; i < nrA; ++i)
	{
		startAi = i*ncA;
		for (j=0; j < ncA; ++j, ++startAi)
		{
			startBi = i*ncB; startCj = j;
			_domain.assign(tmp, _domain.zero);
			for (k=0; k<ncB; ++k, ++startBi, startCj+=ncA) // ncA = ncC
				_domain.axpy(tmp, B[startBi], C[startCj], tmp);
			_domain.assign( A[startAi], tmp );
		}
	}
}

template<class Domain>
void MatrixDom<Domain,Dense>::mulin
( Rep& R, const Type_t& val) const
{
	_supportdomain.mulin(R, val);
}

template<class Domain>
void MatrixDom<Domain,Dense>::mul
( typename VectorDom<Domain,Dense>::Rep& res,
  const Rep& M,
  const VectorDom<Domain,Dense>& VD,
  const typename VectorDom<Domain,Dense>::Rep& u ) const
{
	GIVARO_ASSERT( (nrow(M) == VD.dim(res)), "Bad size");
	GIVARO_ASSERT( (ncol(M) == VD.dim(u)), "Bad size");
	Indice_t k, i, j;
	int startMi;
	Indice_t nrM = nrow(M);
	Indice_t ncM = ncol(M);

	// -- the algorithm used the fact that matrices are stored by row
	for (i=0; i < nrM; ++i)
	{
		_domain.assign(res[i], _domain.zero);
		startMi = i*ncM;
		for (j=0; j < ncM; ++j, ++startMi)
			_domain.axpyin(res[i], M[startMi], u[j]);
	}
}


// res <- tr(u) * tr(A)
template<class Domain>
void MatrixDom<Domain,Dense>::multrans
( typename VectorDom<Domain,Dense>::Rep& res,
  const Rep& M,
  const VectorDom<Domain,Dense>& VD,
  const typename VectorDom<Domain,Dense>::Rep& u ) const
{
	GIVARO_ASSERT( (nrow(M) == VD.dim(u)), "Bad size");
	GIVARO_ASSERT( (ncol(M) == VD.dim(res)), "Bad size");
	Indice_t k, i, j;
	int startMi;
	Indice_t nrM = nrow(M);
	Indice_t ncM = ncol(M);

	// -- the algorithm used the fact that matrices are stored by row
	for (j=0; j < ncM; ++j)
		_domain.assign(res[j], _domain.zero);
	for (i=0; i < nrM; ++i)
	{
		startMi = i*ncM;
		for (j=0; j < ncM; ++j, ++startMi)
			_domain.axpyin(res[j], u[i], M[startMi]);
	}
}



template<class Domain>
void MatrixDom<Domain,Dense>::mul
( Rep& R, const Type_t& val, const Rep& A) const
{
	GIVARO_ASSERT((ncol(R) == ncol(A)), "Bad size");
	GIVARO_ASSERT((nrow(R) == nrow(A)), "Bad size");
	_supportdomain.mul(R, val, A);
}

template<class Domain>
void MatrixDom<Domain,Dense>::mul
( Rep& R, const Rep& A, const Type_t& val) const
{
	GIVARO_ASSERT((ncol(R) == ncol(A)), "Bad size");
	GIVARO_ASSERT((nrow(R) == nrow(A)), "Bad size");
	_supportdomain.mul(R, A, val);
}

// r <- a*x+y
template<class Domain>
void MatrixDom<Domain,Dense>::axpy
( Rep& R, const Type_t& a, const Rep& X, const Rep& Y ) const
{
	GIVARO_ASSERT((ncol(R) == ncol(X)) && (ncol(X) == ncol(Y)), "Bad size");
	GIVARO_ASSERT((nrow(R) == nrow(X)) && (nrow(X) == nrow(Y)), "Bad size");
	_supportdomain.axpy(R, a, X, Y);
}
// r <- r+a*x
template<class Domain>
void MatrixDom<Domain,Dense>::axpyin
( Rep& R, const Type_t& a, const Rep& X ) const
{
	GIVARO_ASSERT((ncol(R) == ncol(X)), "Bad size");
	GIVARO_ASSERT((nrow(R) == nrow(X)), "Bad size");
	_supportdomain.axpyin(R, a, X);
}
// r <- a*x-y
template<class Domain>
void MatrixDom<Domain,Dense>::axmy ( Rep& R,
				     const Type_t& a, const Rep& X, const Rep& Y ) const
{
	GIVARO_ASSERT((ncol(R) == ncol(X)) && (ncol(X) == ncol(Y)), "Bad size");
	GIVARO_ASSERT((nrow(R) == nrow(X)) && (nrow(X) == nrow(Y)), "Bad size");
	_supportdomain.axmy(R, a, X, Y);
}
// r <- a*x-r
template<class Domain>
void MatrixDom<Domain,Dense>::axmyin ( Rep& R,
				       const Type_t& a, const Type_t& X ) const
{
	GIVARO_ASSERT((ncol(R) == ncol(X)), "Bad size");
	GIVARO_ASSERT((nrow(R) == nrow(X)), "Bad size");
	_supportdomain.axmy(R, a, X);
}


// -- axpy operations:
// A*X + Y
template<class Domain>
void MatrixDom<Domain,Dense>::axpy
( Rep& R, const Rep& A, const Rep& X, const Rep& Y ) const
{
	int k, i, j;
	int startAi, startXj, startRi;
	int nrR = nrow(R);
	int ncR = ncol(R);
	int ncA = ncol(A);

	// -- the algorithm used the fact that matrices are stored by row
	Type_t tmp;
	for (i=0; i < nrR; ++i)
	{
		startRi = i*ncR;
		for (j=0; j < ncR; ++j, ++startRi)
		{
			startAi = i*ncA; startXj = j;
			_domain.assign(tmp, _domain.zero);
			for (k=0; k<ncA; ++k, ++startAi, startXj+=ncR) // ncR = ncX
				_domain.axpy(tmp, A[startAi], X[startXj], tmp);
			_domain.add( R[startRi], tmp, Y[startRi] );
		}
	}
}

// a*A*X - b*Y
template<class Domain>
void MatrixDom<Domain,Dense>::axpy
( Rep& res, const Type_t& a, const Rep& A,
  const Rep& X, const Type_t& b, const Rep& Y ) const
{
	Indice_t k, i, j;
	int startAi, startXj, startRi;
	Indice_t nrR = nrow(R);
	Indice_t ncR = ncol(R);
	Indice_t ncA = ncol(A);

	// -- the algorithm used the fact that matrices are stored by row
	Type_t tmp;
	for (i=0; i < nrR; ++i)
	{
		startRi = i*ncR;
		for (j=0; j < ncR; ++j, ++startRi)
		{
			startAi = i*ncA; startXj = j;
			_domain.assign(tmp, _domain.zero);
			for (k=0; k<ncA; ++k, ++startAi, startXj+=ncR) // ncR = ncX
				_domain.axpy(tmp, A[startAi], X[startXj], tmp);
			_domain.mul( tmp, a, tmp );
			_domain.axpy( R[startRi], b, Y[startRi], tmp);
		}
	}
}

// Map: used for constant operation !, not that OP
template<class Domain>
template<class OP>
void MatrixDom<Domain,Dense>::map
( Rep& R, OP& op ) const
{
	Indice_t i;
	Indice_t max = _supportdomain.dim(R);
	for (i=max; --i>=0; ) op(R[i]);
}

template<class Domain>
template<class OP>
void MatrixDom<Domain,Dense>::map
( Rep& R, OP& op, const Rep& A ) const
{
	Indice_t i;
	Indice_t max = _supportdomain.dim(R);
	for (i=max; --i>=0; ) op(R[i], A[i]);
}

template<class Domain>
template<class OP>
void MatrixDom<Domain,Dense>::map
( Rep& R, OP& op, const Rep& A, const Rep& B ) const
{
	Indice_t i;
	Indice_t max = _supportdomain.dim(R);
	for (i=max; --i>=0; ) op(R[i], A[i], B[i]);
}


template <class Domain>
inline ostream& MatrixDom<Domain,Dense>::write( ostream& sout ) const
{
	return _domain.write(sout << '(') << ",Dense)";
}


// -- read the domain
template<class Domain>
istream& MatrixDom<Domain,Dense>::read( istream& sin )
{
	char ch;
	sin >> std::ws >> ch;
	if (ch != '(')
		GivError::throw_error(
				      GivBadFormat("MatrixDom<Domain,Dense>::read: syntax error no '('"));

	// -- read subdomain:
	_domain.read(sin);

	// -- read ,
	sin >> std::ws >> ch;
	if (ch != ',')
		GivError::throw_error(
				      GivBadFormat("MatrixDom<Domain,Dense>::read: syntax error no ','"));

	// -- read dense:
	char name[6];
	sin >> std::ws; sin.read(name, 5); name[5] = (char)0;
	if (strcmp(name, "Dense") !=0)
		GivError::throw_error(
				      GivBadFormat("MatrixDom<Domain,Dense>::read: syntax error no 'Dense'"));

	// -- read )
	sin >> std::ws >> ch;
	if (ch != ')')
		GivError::throw_error(
				      GivBadFormat("MatrixDom<Domain,Dense>::read: syntax error no ')'"));
	return sin;
}


template <class Domain>
ostream& MatrixDom<Domain,Dense>::write( ostream& sout, const Rep& A) const
{
	int i,j;
	sout << '[' << nrow(A) << ',' << ncol(A) << ",[";
	for (i=0; i < nrow(A); i++)
	{
		sout << "[";
		for (j=0; j < ncol(A); j++)
		{
			_domain.write(sout,A(i,j));
			if (j < ncol(A) - 1) sout << ",";
		}
		if (i < nrow(A) - 1) sout << "]," << endl;
		else sout << "]";
	}
	return sout << "]]";
}

template <class Domain>
istream& MatrixDom<Domain,Dense>::read (istream& sin, Rep& R) const
{
	char ch;
	Indice_t i,j;
	long nr,nc;

	//  -- read [
	sin >> std::ws >> ch;
	if (ch != '[') {
		cerr << "Read: " << ch << endl;
		GivError::throw_error(
				      GivBadFormat("MatrixDom<Domain,Dense>::read: syntax error no '['"));
	}

	// -- read dimensions: nr , nc
	// -- read nrow
	sin >> std::ws >> nr;
	if (nr < 0)
		GivError::throw_error(
				      GivBadFormat("MatrixDom<Domain,Dense>::read: bad row dimension"));

	//  -- read ,
	sin >> std::ws >> ch;
	if (ch != ',')
		GivError::throw_error(
				      GivBadFormat("MatrixDom<Domain,Dense>::read: syntax error no ','"));

	//  -- read ncol
	sin >> std::ws >> nc;
	if (nc < 0)
		GivError::throw_error(
				      GivBadFormat("MatrixDom<Domain,Dense>::read: bad col dimension"));

	//  -- read ,
	sin >> std::ws >> ch;
	if (ch != ',')
		GivError::throw_error(
				      GivBadFormat("MatrixDom<Domain,Dense>::read: syntax error no ','"));

	//  -- read [
	sin >> std::ws >> ch;
	if (ch != '[')
		GivError::throw_error(
				      GivBadFormat("MatrixDom<Domain,Dense>::read: syntax error no '[' before entries"));

	// -- Read the matrix:
	R.resize( nr, nc );
	for (i=0; i<nr; ++i) {
		if (i != 0) {
			sin >> std::ws >> ch;
			if (ch != ',')
				GivError::throw_error(
						      GivBadFormat("MatrixDom<Domain,Dense>::read: syntax error no ','"));
		}

		cout << "Row: " << i << endl;
		//  -- read [
		sin >> std::ws >> ch;
		if (ch != '[')
			GivError::throw_error(
					      GivBadFormat("MatrixDom<Domain,Dense>::read: syntax error no '['"));

		// read nc-1 times: val ,
		for (j=0; j< nc-1; ++j) {
			cout << "Read: i=" << i << ", j=" << j << endl;
			_domain.read(sin, R(i,j));
			//  -- read [
			sin >> std::ws >> ch;
			if (ch != ',')
				GivError::throw_error(
						      GivBadFormat("MatrixDom<Domain,Dense>::read: syntax error no ','"));
		}

		// read 1 : val ]
		cout << "Read: i=" << i << ", j=" << j << endl;
		_domain.read(sin, R(i,j));
		sin >> std::ws >> ch;
		if (ch != ']')
			GivError::throw_error(
					      GivBadFormat("MatrixDom<Domain,Dense>::read: syntax error no ']'"));
	}

	//  -- read ] ]
	sin >> std::ws >> ch;
	if (ch != ']')
		GivError::throw_error(
				      GivBadFormat("MatrixDom<Domain,Dense>::read: syntax error no ']'"));
	sin >> std::ws >> ch;
	if (ch != ']')
		GivError::throw_error(
				      GivBadFormat("MatrixDom<Domain,Dense>::read: syntax error no ']'"));
	return sin;
}
} // Givaro

//#include "givaro/givmatdenseops.f.spe"
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
