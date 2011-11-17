// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givpoly1muldiv.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givpoly1muldiv.inl,v 1.14 2011-02-02 16:23:56 bboyer Exp $
// ==========================================================================
#ifndef __GIVARO_poly1_muldiv_INL
#define __GIVARO_poly1_muldiv_INL
#include "givaro/givpower.h"
#include "givaro/giverror.h"

namespace Givaro {

template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::shiftin ( Rep& R, int s) const
{
	R.insert(R.begin(), s, this->_domain.zero );
	return R;
}

template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::shift ( Rep& R, const Rep& a, int s) const
{
	R = a;
	return R.shiftin(R, s);
}


template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::mulin( Rep& R, const Type_t& u ) const
{
	for(typename Rep::iterator ri = R.begin();ri!=R.end();++ri)
		_domain.mulin(*ri, u);
	return R;

	//  return _supportdomain.mulin(R,u);
}

template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::mulin( Rep& R, const Rep& P ) const
{
	size_t sR = R.size();
	size_t sP = P.size();
	Rep tmp(sR+sP);
	mul(tmp, R, P);
	//   R.logcopy(tmp);
	//   return R;
	return assign(R,tmp);
}

#if 0
template <class Domain>
Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::mul( Rep& R, const Rep& P, const Rep& Q ) const
{
	size_t sR = R.size();
	size_t sP = P.size();
	size_t sQ = Q.size();
	if ((sQ ==0) || (sP ==0)) { R.reallocate(0); return R; }
	if (sR != sQ+sP) R.reallocate(sR = sP+sQ-1);

	size_t i,j;
	for (i=0; i<sR; ++i) _domain.assign(R[i], _domain.zero);
	for (i=0; i<sP; ++i)
		if (! _domain.isZero(P[i]))
			for (j=0; j<sQ; ++j)
				_domain.axpy(R[i+j], P[i], Q[j], R[i+j]);
	return setdegree(R);
}
#endif

template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::mul( Rep& R, const Rep& P, const Rep& Q ) const
{
	size_t sR = R.size();
	size_t sP = P.size();
	size_t sQ = Q.size();
	if ((sQ ==0) || (sP ==0)) { R.reallocate(0); return R; }
	if (sR != sQ+sP) R.reallocate(sR = sP+sQ-1);

//  std::cerr << "MULBASE rbeg: " << &( *(R.begin())) << std::endl;
//  std::cerr << "MULBASE rend: " << (&(R.back())+1) << std::endl;


 	mul(R, R.begin(), R.end(),
            P, P.begin(), P.end(),
            Q, Q.begin(), Q.end());
        
        return setdegree(R);
}

template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::stdmul( Rep& R, const Rep& P, const Rep& Q ) const
{
	size_t sR = R.size();
	size_t sP = P.size();
	size_t sQ = Q.size();
	if ((sQ ==0) || (sP ==0)) { R.reallocate(0); return R; }
	if (sR != sQ+sP) R.reallocate(sR = sP+sQ-1);

//  std::cerr << "MULBASE rbeg: " << &( *(R.begin())) << std::endl;
//  std::cerr << "MULBASE rend: " << (&(R.back())+1) << std::endl;


 	stdmul(R, R.begin(), R.end(),
            P, P.begin(), P.end(),
            Q, Q.begin(), Q.end());
        
        return setdegree(R);
}

#ifndef KARA_THRESHOLD
#define KARA_THRESHOLD 50
#endif

template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::mul( Rep& R, const typename Rep::iterator Rbeg, const typename Rep::iterator Rend, const Rep& P, const typename Rep::const_iterator Pbeg, const typename Rep::const_iterator Pend, const Rep& Q, const typename Rep::const_iterator Qbeg, const typename Rep::const_iterator Qend ) const
{
//  std::cerr << "MULptr rbeg: " << &( *(Rbeg)) << std::endl;
//  std::cerr << "MULptr rend: " << &(*(Rend)) << std::endl;
//  std::cerr << "MULptr r.beg: " << &( *(R.begin())) << std::endl;
//  std::cerr << "MULptr r.end: " << (&(R.back())+1) << std::endl;

    if ( ( (Pend-Pbeg)> KARA_THRESHOLD ) &&
         ( (Qend-Qbeg)> KARA_THRESHOLD) )
        return karamul(R, Rbeg, Rend,
                       P, Pbeg, Pend,
                       Q, Qbeg, Qend);
    else 
        return stdmul(R, Rbeg, Rend,
                      P, Pbeg, Pend,
                      Q, Qbeg, Qend);
}

template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::stdmul( Rep& R, const typename Rep::iterator Rbeg, const typename Rep::iterator Rend, const Rep& P, const typename Rep::const_iterator Pbeg, const typename Rep::const_iterator Pend, const Rep& Q, const typename Rep::const_iterator Qbeg, const typename Rep::const_iterator Qend ) const
{

// std::cerr << "-----------------------" << std::endl;
// //  std::cerr << "MULkar rbeg: " << &( *(Rbeg)) << std::endl;
// //  std::cerr << "MULkar rend: " << &(*(Rend)) << std::endl;
// //  std::cerr << "MULkar r.beg: " << &( *(R.begin())) << std::endl;
// //  std::cerr << "MULkar r.end: " << (&(R.back())+1) << std::endl;

//  this->write(std::cerr << P.size() << '[', P) << "]_" << (Pbeg-P.begin()) << '^' << (Pend-P.begin()) << std::endl;
//  this->write(std::cerr << Q.size() << '[', Q) << "]_" << (Qbeg-Q.begin()) << '^' << (Qend-Q.begin()) << std::endl;
//  affiche(std::cerr << '[', _domain, R) << "]^" << (Rend-Rbeg) << std::endl;
    
// //  std::cerr << "MULkar2 rbeg: " << &( *(Rbeg)) << std::endl;
// //  std::cerr << "MULkar2 rend: " << &(*(Rend)) << std::endl;
// //  std::cerr << "MULkar2 r.beg: " << &( *(R.begin())) << std::endl;
// //  std::cerr << "MULkar2 r.end: " << (&(R.back())+1) << std::endl;

	typename Rep::const_iterator ai=Pbeg,bi=Qbeg;
	typename Rep::iterator ri=Rbeg, rig=Rbeg;
	if (_domain.isZero(*ai))
		for(;bi!=Qend;++bi,++ri)
// { 
			*ri = _domain.zero;
// _domain.write(std::cerr << "r00: ", *ri) << std::endl;
// affiche(std::cerr << "R: ", _domain, R) << std::endl;
// }
        

	else
		for(;bi!=Qend;++bi,++ri)
// { 
			if (_domain.isZero(*bi))
				*ri = _domain.zero;
			else
				_domain.mul(*ri,*ai,*bi);
//  std::cerr << "MULkar3 rbeg: " << &( *(Rbeg)) << std::endl;
//  std::cerr << "MULkar3 rend: " << &(*(Rend)) << std::endl;
//  std::cerr << "MULkar3 r.beg: " << &( *(R.begin())) << std::endl;
//  std::cerr << "MULkar3 r.end: " << (&(R.back())+1) << std::endl;

// _domain.write(std::cerr << "r0a: ", *ri) << std::endl;
//  std::cerr << "rbeg: " << &( *(R.begin())) << std::endl;
//  std::cerr << "ri: " << &(*ri) << std::endl;
//  std::cerr << "rend: " << (&(R.back())+1) << std::endl;

// affiche(std::cerr << "R: ", _domain, R) << std::endl;
// }
        
	for(;ri!=Rend;++ri)
		*ri = _domain.zero;
	for(++ai,++rig;ai!=Pend;++ai,++rig)
		if (! _domain.isZero(*ai))
			for(ri=rig,bi=Qbeg;bi!=Qend;++bi,++ri)
// { 
				_domain.axpyin(*ri,*ai,*bi);
// _domain.write(std::cerr << "ri: ", *ri) << std::endl;
// affiche(std::cerr << "R: ", _domain, R) << std::endl;
// }
// affiche(std::cerr << "R: ", _domain, R) << std::endl;
//  std::cerr << "MULkar4 rbeg: " << &( *(Rbeg)) << std::endl;
//  std::cerr << "MULkar4 rend: " << &(*(Rend)) << std::endl;
//  std::cerr << "MULkar4 r.beg: " << &( *(R.begin())) << std::endl;
//  std::cerr << "MULkar4 r.end: " << (&(R.back())+1) << std::endl;
// affiche(std::cerr << "R: ", _domain, R) << std::endl;
// std::cerr << "-----------------------" << std::endl;
	return R;
}

#ifndef GIVMIN
#define GIVMIN(a,b) ((a)<(b)?(a):(b))
#endif

template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::karamul( Rep& R, const typename Rep::iterator Rbeg, const typename Rep::iterator Rend, const Rep& P, const typename Rep::const_iterator Pbeg, const typename Rep::const_iterator Pend, const Rep& Q, const typename Rep::const_iterator Qbeg, const typename Rep::const_iterator Qend ) const
{
// std::cerr << "-------------- BEG KARAMUL ----------" << (Pend-Pbeg) << std::endl;
//     std::cerr << "sP: " << (Pend-Pbeg) << std::endl;
//     std::cerr << "sQ: " << (Qend-Qbeg) << std::endl;
//     std::cerr << "sR: " << (Rend-Rbeg) << std::endl;
    for(typename Rep::iterator ri=Rbeg; ri!= Rend; ++ri) _domain.assign(*ri,_domain.zero);
    

    size_t halfP = (Pend-Pbeg)>>1;
    size_t halfQ = (Qend-Qbeg)>>1;
    size_t half = GIVMIN(halfP, halfQ);
    size_t halfR = half<<1;
//     std::cerr << "hP: " << halfP << std::endl;
//     std::cerr << "hQ: " << halfQ << std::endl;
//     std::cerr << "h1R: " << half << std::endl;
//     std::cerr << "hR: " << halfR << std::endl;

    typename Rep::const_iterator Pmid=Pbeg+half;
    typename Rep::const_iterator Qmid=Qbeg+half;
    typename Rep::iterator Rmid=Rbeg+halfR;
    
//     this->write(std::cerr << "P:=", P) << ';' << std::endl;
//     this->write(std::cerr << "Q:=", Q) << ';' << std::endl;

    mul(R, Rbeg, Rmid, 
        P, Pbeg, Pmid,
        Q, Qbeg, Qmid);
    
//     affiche(std::cerr << "PlQl:=", _domain, R) << ';' << std::endl;

    mul(R, Rmid, Rend,
        P, Pmid, Pend,
        Q, Qmid, Qend);
    
//     affiche(std::cerr << "PhQh:=", _domain, R) << ';' << std::endl;

    Rep PHPL;
    for(typename Rep::const_iterator PHi=Pmid; PHi!=Pend; ++PHi)
        PHPL.push_back(*PHi);
//     affiche(std::cerr << "Ph:=", _domain, PHPL) << ';' << std::endl;
    subin(PHPL, PHPL.begin(), P, Pbeg, Pmid); // Ph - Pl
    setdegree(PHPL);
//     affiche(std::cerr << "Ph-Pl: ", _domain, PHPL) << std::endl;
    
    Rep QHQL;
    for(typename Rep::const_iterator QHi=Qmid; QHi!=Qend; ++QHi)
        QHQL.push_back(*QHi);
//     affiche(std::cerr << "Qh:=", _domain, QHQL) << ';' << std::endl;
    subin(QHQL, QHQL.begin(), Q, Qbeg, Qmid); // Qh - Ql
    setdegree(QHQL);
//     affiche(std::cerr << "Qh-Ql: ", _domain, QHQL) << std::endl;
    
    Rep M;
    mul(M, PHPL, QHQL);
    setdegree(M);
    
//     affiche(std::cerr << "M1:=", _domain, M) << ';' << std::endl;
//     std::cerr << "Rmid-Rbeg: " << (Rmid-Rbeg) << std::endl;
    subin(M, M.begin(), M.end(), R, Rbeg, Rmid); // -= PlQl
    setdegree(M);
//     affiche(std::cerr << "M2:=", _domain, M) << ';' << std::endl;
    subin(M, M.begin(), M.end(), R, Rmid, Rend);		// -= PhQh
    setdegree(M);
//     affiche(std::cerr << "M3:=", _domain, M) << ';' << std::endl;


    typename Rep::iterator ri=Rbeg+half;
    typename Rep::const_iterator mi=M.begin();
    for( ; mi != M.end(); ++ri, ++mi) _domain.subin(*ri, *mi);

//     affiche(std::cerr << "R:=", _domain, R) << ';' << std::endl;
//     std::cerr << "-------------- END KARAMUL ----------" << (Pend-Pbeg) << std::endl;
    return R;
}




// Compute truncated mul: only the coefficients inside the degree interval, included
template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::mul( Rep& R, const Rep& P, const Rep& Q, const Degree& Val, const Degree& deg) const
{
	size_t sR = R.size();
	size_t sP = P.size();
	size_t sQ = Q.size();
	if ((sQ ==0) || (sP ==0)) { R.reallocate(0); return R; }
	size_t newS = value(deg-Val)+1;
	if (sR != newS) R.reallocate(sR = newS);
	for(typename Rep::iterator ri=R.begin(); ri!= R.end(); ++ri)
		*ri = _domain.zero;

	for(size_t i=0; i<sR; ++i) {
		long k=i+Val.value();
		size_t j=0;
		if (static_cast<size_t>(k)>=sQ) {
			j=k;
			k=sQ-1;
			j-=k;
		}
		for( ; (j<sP) && (k>=0); ++j,--k) {
			_domain.axpyin(R[i],P[j],Q[k]);
		}
	}
	return setdegree(R);
}

template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::mul
( Rep& R, const Rep& P, const Type_t& u ) const
{
	typename Rep::const_iterator ip = P.begin();
	R.resize(P.size());
	for(typename Rep::iterator ir = R.begin(); ir != R.end(); ++ir, ++ip)
		this->_domain.mul(*ir, *ip, u);
	return R;
}

template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::mul
( Rep& R, const Type_t& u, const Rep& P ) const
{
	return this->mul(R,P,u);
}
} // Givaro

//#include <typeinfo>

namespace Givaro {

template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::divin(Rep& R, const Type_t& u) const
{
#ifdef GIVARO_DEBUG
	if (_domain.isZero(u)) GivError::throw_error(GivMathDivZero("[Poly1Dom<D>::divin]"));
#endif
	size_t sz =R.size();
	for (unsigned int i=0; i<sz; ++i)
		_domain.divin(R[i],u);
	return setdegree(R);
}


template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::div(Rep& R, const Rep& P, const Type_t& u) const
{
#ifdef GIVARO_DEBUG
	if (_domain.isZero(u)) GivError::throw_error(GivMathDivZero("[Poly1Dom<D>::div]"));
#endif
	size_t sP =P.size();
	R.reallocate(sP);
	for (unsigned int i=0; i<sP; ++i)
		_domain.div(R[i],P[i],u);
	return setdegree(R);
}


template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::div(Rep& R, const Type_t& u, const Rep& P) const
{
#ifdef GIVARO_DEBUG
	if (isZero(P)) GivError::throw_error(GivMathDivZero("[Poly1Dom<D>::div]"));
#endif
	if (_domain.isZero(u)) { return assign(R,zero);}
	size_t sP =P.size();
	if (sP >1) { R.reallocate(0); return R; }
	size_t sR =R.size();
	if (sR !=1) R.reallocate(1);
	_domain.div(R[0], u, P[0]);
	return setdegree(R);
}

template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::div(Rep& Q, const Rep& A, const Rep& B) const
{
	Rep R;
	return divmod(Q,R,A,B);
}

template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::divin(Rep& Q, const Rep& A) const
{
	Rep R, B;
	divmod(B,R,Q,A);
	return assign(Q,B);
}

template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::modin(Rep& R, const Type_t& u) const
{
#ifdef GIVARO_DEBUG
	if (_domain.isZero(u)) GivError::throw_error(GivMathDivZero("[Poly1Dom<D>::modin]"));
#endif
	R.reallocate(0);
	return R;
}


template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::mod(Rep& R, const Rep& P, const Type_t& u) const
{
#ifdef GIVARO_DEBUG
	if (_domain.isZero(u)) GivError::throw_error(GivMathDivZero("[Poly1Dom<D>::mod]"));
#endif
	R.reallocate(0);
	return R;
}


template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::mod(Rep& R, const Type_t& u, const Rep& P) const
{
#ifdef GIVARO_DEBUG
	if (isZero(P)) GivError::throw_error(GivMathDivZero("[Poly1Dom<D>::mod]"));
#endif
	if (_domain.isZero(u)) { return assign(P,R); }
	size_t sP =P.size();
	if (sP >1) {
		R.reallocate(1);
		_domain.assign(R[0], u);
		return R;
	}
	R.reallocate(1);
	_domain.mod(R[0],u,P[0]);
	return R;
}

#if 0
template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::modin(Rep& R, const Rep& A) const
{
	Rep tR; assign(tR,R);
	return mod(R,tR,A);
}
#endif

template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::modin(Rep& A, const Rep& B) const
{
	// In place remainder
	// A is written with next remainder in
	// the division algorithm written at the end.
	// Last step is erasing of the first values.
	//     write(std::cerr << "Rem(", A) << " ,";
	//     write(std::cerr, B) << ", X) mod " << _domain.size();
	long i = A.size()-B.size();
	if (i >= 0) {
		typedef typename Rep::value_type TT;
		TT l;
		typename Rep::reverse_iterator ai,aai;
		typename Rep::const_reverse_iterator bi;
		for (; i>=0; --i) {
			ai = A.rbegin();
			bi = B.rbegin();
			_domain.div(l,*ai,*bi);
			aai = A.rbegin();
			for(++bi,++ai;bi!=B.rend();++bi,++ai,--i) {
				_domain.maxpy(*aai,l,*bi,*ai);
				if (! _domain.isZero(*aai)) break;
			}
			if (bi!=B.rend())
				for(++bi,++ai,++aai;bi!=B.rend();++bi,++ai,++aai)
					_domain.maxpy(*aai,l,*bi,*ai);
			for(;ai!=A.rend();++ai,++aai)
				*aai = *ai;
			*aai = _domain.zero;
		}
		//         write(std::cerr << " = ", A) << ";" << std::endl;
		A.erase(A.begin(), A.begin()+(A.size()-B.size()-i));
	}
	//     write(std::cerr << " = ", setdegree(A)) << ";" << std::endl;
	return setdegree(A);
}

template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::mod(Rep& R, const Rep& A, const Rep& B) const
{
	Rep Q;
	//   write(std::cerr, A) << " = (";
	//   write(std::cerr, B) << ") * (";
	divmod(Q,R,A,B);
	//   write(std::cerr, Q) << ") + (";
	//   write(std::cerr, R) << ");" << std::endl;
	return R;
}

// #include <typeinfo>

template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::divmod( Rep& Q, Rep& R, const Rep& A,  const Rep& B) const
// returns Q such that A = B Q + R
{
	//     std::cerr << "BEG divmod of " << typeid(*this).name() << std::endl;
	//     std::cerr << "BEG with _domain " << typeid(_domain).name() << std::endl;
	Degree degB;
	degree(degB, B);
#ifdef GIVARO_DEBUG
	if (degB == Degree::deginfty)
		GivError::throw_error(GivMathDivZero("[Poly1Dom<D>::div]"));
#endif
	Degree degA; degree(degA, A);
	if (degA == Degree::deginfty) {
		assign(R, zero);
		return assign(Q, zero);
	}
	if (degB == 0) // cste
	{
		assign(R, zero);
		return div(Q, A, B[0]);
	}
	//   write(std::cerr, A) << " of degA " << degA << std::endl;
	//   write(std::cerr, B) << " of degB " << degB << std::endl;


	// JGD 15.12.1999 :
	//   if (degA ==0)
	//   {
	//     assign(R, zero);
	//     return assign(Q, zero);
	//   }
	if (degB > degA) {
		assign(R, A);
		return assign(Q, zero);
	}

	long degQuo = value(degA-degB);
	long degRem = value(degA);
	Q.reallocate(degQuo+1);

	assign(R,A);
	// write(std::cerr << "A:=", A) << "; # of degA " << degA << std::endl;
	// write(std::cerr << "B:=", B) << "; # of degB " << degB << std::endl;

	for (long i=degQuo; i>=0; --i)
	{
		// == ld X^ (degRem-degQ)
		_domain.div(Q[i], R[degRem], B[degB.value()]);
		// _domain.write(std::cerr << "Q[" << i << "]:=", Q[i]) << ';' << std::endl;
		//  std::cerr << "degB: " << degB << std::endl;
		for (long j=0; degB>j; ++j) { // rem <- rem - ld*x^(degRem-degB)*B
			_domain.maxpyin(R[j+i], Q[i], B[j]);
		}
		_domain.assign(R[degRem],_domain.zero) ; --degRem;
		// write(std::cerr << "inR:=", R) << ';' << std::endl;
	}
	// write(std::cerr << "Q:=", Q) << "; # of degQ " << degQuo << std::endl;
	// write(std::cerr << "R:=", R) << "; # of degR " << degRem << std::endl;
	R.reallocate(degRem+1);
	setdegree(R);
	//     std::cerr << "END divmod of " << typeid(*this).name() << std::endl;
	return setdegree(Q);
}

template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::divmodin( Rep& Q, Rep& R, const Rep& B) const
// returns Q such that A = B Q + R
{
	//     std::cerr << "BEG divmod of " << typeid(*this).name() << std::endl;
	//     std::cerr << "BEG with _domain " << typeid(_domain).name() << std::endl;
	Degree degB; degree(degB, B);
#ifdef GIVARO_DEBUG
	if (degB == Degree::deginfty)
		GivError::throw_error(GivMathDivZero("[Poly1Dom<D>::div]"));
#endif
	Degree degA; degree(degA, R);
	if (degA == Degree::deginfty) {
		assign(R, zero);
		return assign(Q, zero);
	}
	if (degB == 0) // cste
	{
		div(Q, R, B[0]);
 		assign(R, zero);
                return Q;
	}
	//   write(std::cerr, A) << " of degA " << degA << std::endl;
	//   write(std::cerr, B) << " of degB " << degB << std::endl;


	// JGD 15.12.1999 :
	//   if (degA ==0)
	//   {
	//     assign(R, zero);
	//     return assign(Q, zero);
	//   }
	if (degB > degA) {
		return assign(Q, zero);
	}

	long degQuo = value(degA-degB);
	long degRem = value(degA);
	Q.reallocate(degQuo+1);

	// write(std::cerr << "A:=", A) << "; # of degA " << degA << std::endl;
	// write(std::cerr << "B:=", B) << "; # of degB " << degB << std::endl;

	for (long i=degQuo; i>=0; --i)
	{
		// == ld X^ (degRem-degQ)
		_domain.div(Q[i], R[degRem], B[degB.value()]);
		// _domain.write(std::cerr << "Q[" << i << "]:=", Q[i]) << ';' << std::endl;
		//  std::cerr << "degB: " << degB << std::endl;
		for (long j=0; degB>j; ++j) { // rem <- rem - ld*x^(degRem-degB)*B
			_domain.maxpyin(R[j+i], Q[i], B[j]);
		}
		_domain.assign(R[degRem],_domain.zero) ; --degRem;
		// write(std::cerr << "inR:=", R) << ';' << std::endl;
	}
	// write(std::cerr << "Q:=", Q) << "; # of degQ " << degQuo << std::endl;
	// write(std::cerr << "R:=", R) << "; # of degR " << degRem << std::endl;
	R.reallocate(degRem+1);
	setdegree(R);
	//     std::cerr << "END divmod of " << typeid(*this).name() << std::endl;
	return setdegree(Q);
}


template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::pdivmod
( Rep& Q, Rep& R, Type_t& m, const Rep& A, const Rep& B) const
// returns Q ...
{
	Degree degB; degree(degB, B);
#ifdef GIVARO_DEBUG
	if (degB == Degree::deginfty)
		GivError::throw_error(GivMathDivZero("[Poly1Dom<D>::div]"));
#endif
	Degree degA; degree(degA, A);
	if (degA == Degree::deginfty) {
		assign(R, zero);
		_domain.assign(m, _domain.one);
		return assign(Q, zero);
	}
	if (degB == 0) // cste
	{
		assign(R, zero);
		_domain.assign(m, B[0]);
		return assign(Q, A);
	}
	if (degA ==0)
	{
		assign(R, zero);
		_domain.assign(m, _domain.one);
		return assign(Q, zero);
	}
	if (degB > degA) {
		assign(R, A);
		_domain.assign(m, _domain.one);
		return assign(Q, zero);
	}

	long degQuo = value(degA-degB);
	long degRem = value(degA);
	Q.reallocate(degQuo+1);
	assign(R,A);

	Type_t tmp, lB;
	_domain.assign(lB, B[degB.value()]);
	_domain.assign(m, _domain.one);
	long i,j;
	for (i=degQuo; i>=0; --i)
	{
		// == ld X^ (degRem-degQ)
		_domain.assign(Q[degQuo], R[degRem]);

		// rem <- lB*rem - lQ*x^(degRem-degB)*B
		for (j=0; j<degQuo; j++)
			_domain.mulin (R[j], lB);
		for (j=0; degB>j; j++)
		{
			_domain.mulin(R[j+degQuo], lB);
			_domain.maxpyin(R[j+degQuo], Q[degQuo], B[j]);
		}
		_domain.assign(R[degRem],_domain.zero); degQuo--; degRem--;
		_domain.mulin(m, lB);
	}
	R.reallocate(degRem+1);
	setdegree(R);
	return setdegree(Q);
	//  Poly1Dom<Domain,Dense>::Rep U,V;
	//  assign(U,A);
	//  mulin(U,m);
	//  write(std::cout << "m*A:", U) << std::endl;
	//  mul(U,Q,B);
	//  write(std::cout << "Q*B:", U) << std::endl;
}



template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::pmod
( Rep& R, Type_t& m, const Rep& A, const Rep& B) const
{
	Degree degB; degree(degB, B);
#ifdef GIVARO_DEBUG
	if (degB == Degree::deginfty)
		GivError::throw_error(GivMathDivZero("[Poly1Dom<D>::div]"));
#endif
	Degree degA; degree(degA, A);
	if (degA == Degree::deginfty) {
		_domain.assign(m, _domain.one);
		return assign(R, zero);
	}
	if (degB == 0) // cste
	{
		_domain.assign(m, B[0]);
		return assign(R, zero);
	}
	if (degA ==0)
	{
		_domain.assign(m, _domain.one);
		return assign(R, zero);
	}
	if (degB > degA) {
		_domain.assign(m, _domain.one);
		return assign(R, A);
	}

	Degree degR = degA;
	assign(R,A);

	Type_t tmp, lB;
	_domain.assign(lB, B[degB.value()]);
	//write(std::cout << "B:", B) << std::endl;
	//_domain.write(std::cout << "lB:", lB) << "^" << degA-degB+1 << std::endl;
	//   _domain.pow(m, lB, degA.value()-degB.value()+1);
	dom_power(m, lB, degA.value()-degB.value()+1,_domain);
	//_domain.write(std::cout << "m:", m) << std::endl;
	for (; degB<= degR; )
	{
		long d = degR.value()-degB.value();
		// R <- lB*R - lR*x^(degR-degB)*B
		for (long j=0; degB>j; j++)
		{
			_domain.mulin (R[j+d], lB);
			_domain.maxpyin(R[j+d], R[degR.value()], B[j]);
		}
		for (long j=0; j<d; ++j)
			_domain.mulin (R[j], lB);
		_domain.assign(R[degR.value()],_domain.zero);
		degree(degR, R);
	}
	R.reallocate(degR.value()+1);
	return setdegree(R);
}

} // Givaro
#endif // __GIVARO_poly1_muldiv_INL
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
