// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givpoly1misc.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givpoly1misc.inl,v 1.22 2011-02-02 16:23:56 bboyer Exp $
// ==========================================================================
// Description:

#include <algorithm>

#ifndef __GIVARO_poly_misc_INL
#define __GIVARO_poly_misc_INL

namespace Givaro {
	template<class Domain>
	inline int Poly1Dom<Domain,Dense>::isZero (const Rep& P) const
	{
		setdegree(const_cast<Rep&>(P));
		if (P.size() ==0) return 1;
		if (P.size() ==1) return _domain.isZero(P[0]);
		else return 0;
	}

	template<class Domain>
	inline int Poly1Dom<Domain,Dense>::isOne ( const Rep& P ) const
	{
		setdegree(const_cast<Rep&>(P));
		if (P.size() ==1)
			return _domain.isOne(P[0]);
		else
			return 0;
	}

	template<class Domain>
	inline int Poly1Dom<Domain,Dense>::areEqual (const Rep& P, const Rep& Q) const
	{
		setdegree(const_cast<Rep&>(P));
		setdegree(const_cast<Rep&>(Q));
		if (P.size() != Q.size()) return 0;
		for( typename Element::const_iterator pit = P.begin(), qit = Q.begin();
		     pit != P.end();
		     ++pit, ++qit)
			if ( !_domain.areEqual(*pit, *qit) ) return 0;
		return 1;
	}


	template<class Domain>
	inline int Poly1Dom<Domain,Dense>::areNEqual (const Rep& P, const Rep& Q) const
	{
		setdegree(const_cast<Rep&>(P));
		setdegree(const_cast<Rep&>(Q));
		if (P.size() != Q.size()) return 1;
		for( typename Element::const_iterator pit = P.begin(), qit = Q.begin();
		     pit != P.end();
		     ++pit, ++qit)
			if ( !_domain.areEqual(*pit, *qit) ) return 1;
		return 0;
	}


	// -- Compute the degree of P
	template <class Domain>
	inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::setdegree ( Rep& P ) const
	{
		int sz = (int) (P.size() - 1);
		if (P.size() <= 0) return P.reallocate(0);
		if (_domain.isZero(P[sz]) ==0) {
			return P;
		}
		for (int j=sz; j--; )
			if (_domain.isZero(P[j]) ==0) {
				P.reallocate(j+1);
				return P;
			}
		return P.reallocate(0);
	}

	template <class Domain>
	inline Degree& Poly1Dom<Domain,Dense>::degree(Degree& deg, const Rep& P) const
	{
		int sz = (int) P.size();
		if (sz ==0) {
			return deg = Degree::deginfty;
		}
		if (_domain.isZero(P[sz-1])) {
			setdegree(const_cast<Rep&>(P));
			sz = (int) P.size();
		}
		return deg = (Degree) (sz-1);
	}


	template<class Domain>
	inline typename Poly1Dom<Domain,Dense>::Type_t& Poly1Dom<Domain,Dense>::leadcoef (Type_t& c, const Rep& P) const
	{
		Degree dP;
		degree(dP, P);
		if (dP == Degree::deginfty) return _domain.assign(c, _domain.zero);
		else return _domain.assign(c, P[dP.value()]);
	}


	// -- Returns the i-th coefficients
	template<class Domain>
	inline typename Poly1Dom<Domain,Dense>::Type_t& Poly1Dom<Domain,Dense>::getEntry (Type_t& c, const Degree& i, const Rep& P) const
	{
		Degree dP;
		degree(dP, P);
		if (dP < i) return _domain.assign(c, _domain.zero);
		else return _domain.assign(c, P[i.value()]);
	}




	template <class Domain>
	inline Degree& Poly1Dom<Domain,Dense>::val(Degree& d, const Rep& P) const
	{
		size_t sz = P.size();
		if (sz ==0) {
			return d = Degree::deginfty;
		}
		if (!_domain.isZero(P[0])) {
			return d = 0;
		}
		for (size_t i=1; i<sz; ++i)
		{
			if (!_domain.isZero(P[i])) {
				return d = i;
			}
		}
		return d=0;
	}


	// Horner's scheme for evaluation
	template <class Domain>
	inline typename Poly1Dom<Domain,Dense>::Type_t& Poly1Dom<Domain,Dense>::eval (Type_t& res, const Rep& P, const Type_t& Val) const
	{
		typename Domain::Element tmp;
		_domain.init(tmp,0UL);
		Degree dP ; degree(dP, P);
		if (dP == Degree::deginfty) _domain.assign(res, _domain.zero);
		else {
			_domain.assign(res, P[dP.value()]);
			for (int i = (int)dP.value(); i--; )
				_domain.assign(res,_domain.axpy(tmp, res, Val, P[i]));
		}
		return res;
	}

	template <class Domain>
	inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::diff(Rep& P, const Rep& Q) const
	{
		Degree dQ;
		degree(dQ, Q);
		if ((dQ == Degree::deginfty) || (dQ == 0)) {
			P.reallocate(0);
			return P;
		}
		P.reallocate(dQ.value());
		Type_t cste; _domain.assign(cste, _domain.zero);
		for (int i=0; dQ>i; ++i) {
			_domain.add(cste, cste, _domain.one);
			_domain.mul(P[i], Q[i+1], cste);
		}
		return P;
	}

	template <class Domain>
	inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::pow( Rep& W, const Rep& P, long n) const
	{
		Rep puiss2;
		assign(puiss2, P); // -- P**(2^k)
		Rep tmp;
		assign(W,one);
		unsigned int p;
		p = (unsigned int) (n < 0 ? -n : n); // cannot treats the case of <0 exponent
		while (p != 0) {
			if (p & 0x1) {
				mul( tmp, W, puiss2);
				assign(W,tmp);
			}
			if ((p >>= 1) != 0) {
				mul( tmp, puiss2, puiss2);
				assign(puiss2,tmp);
			}
		}
		return W;
	}


	template <class Domain>
	inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::powmod( Rep& W, const Rep& P, IntegerDom::Element pwr, const Rep& U) const
	{
		IntegerDom ID;
		// ID.write(cerr << "\n----------- POWMOD -----------\n pwr: ", pwr) << endl;
		// write(cerr << "P: ",P) << endl;
		// write(cerr << "U: ",U) << endl;
		Rep puiss, tmp;
		mod(puiss, P, U);
		assign(W,one);
		IntegerDom::Element n,q,r,deux,Zero;
		ID.init(deux,2);
		ID.init(Zero,0UL);
		if (ID.islt(pwr,Zero) )
			ID.neg(n,pwr);
		else
			ID.init(n,pwr);

		while (n > 0) {
			ID.divmod(q,r,n,deux);
			n.copy(q);
			if (! ID.isZero(r)) {
				mulin(W,puiss);
				mod(tmp, W, U);
				assign(W,tmp);
			}
			mul(tmp,puiss,puiss);
			mod(puiss,tmp,U);
		}

		// write(cerr << "W: ", W) << "\n----------- END POWMOD -----------" <<  endl;
		return setdegree(W);
	}


	// -- Random dense polynomial of degree 0
	template <class Domain> template<class RandIter>
	inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::random(RandIter& g, Rep& r) const
	{
		return random(g, r,Degree(0));
	}

	// -- Random dense polynomial of size s
	template <class Domain> template<class RandIter>
	inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::random(RandIter& g, Rep& r, long s) const
	{
		return random(g, r,Degree(s-1));
	}


	// -- Random dense polynomial of degree d
	// TEMPORARY VERSION SATISFYING LINBOX
	template <class Domain> template<class RandIter>
	inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::random(RandIter& g, typename Poly1Dom<Domain,Dense>::Rep& r, Degree d) const
	{
		r.reallocate(d.value()+1);
		typename Domain::Element tmp;
		while (_domain.isZero(g.random(tmp))) ;
		r[d.value()] = tmp;
		for (int i=(int)d.value(); i--;)
			g.random(r[i]);
		return r;
	}

	// -- Random dense polynomial of degree d
	// specialization pour Givrandom
	// #if !defined(__INTEL_COMPILER) && !defined(__CUDACC__)
	// template<>
	// #endif
	template <class Domain>
	inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::random(GivRandom& g, Rep& r, Degree d) const
	{
		r.reallocate(d.value()+1);
		_domain.nonzerorandom(g, r[d.value()]);
		for (int i= (int)d.value(); i--;)
			_domain.random(g,r[i]);
		return r;
	}

#if 0
	template <class Domain> template<class RandIter>
	inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::random(RandIter& g, Rep& r, Degree d) const
	{
		r.reallocate(d.value()+1);
		while (_domain.isZero(_domain.init(r[d.value()], g()))) {};
		for (int i=d.value(); i--;)
			_domain.init(r[i],g());
		_domain.nonzerorandom(g, r[d.value()]);
		for (int i=d.value(); i--;)
			_domain.random(g,r[i]);
		return r;
	}
#endif

	// -- Random dense polynomial with same size as b.
	template <class Domain> template<class RandIter>
	inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::random(RandIter& g, Rep& r, const Rep& b) const
	{
		return random(g, r,b.size());
	}


	template <class Domain> template<class RandIter>
	inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::nonzerorandom(RandIter& g, Rep& r) const
	{
		return random(g, r);
	}


	template <class Domain> template<class RandIter>
	inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::nonzerorandom(RandIter& g, Rep& r, long s) const
	{
		return random(g, r,s);
	}

	template <class Domain> template<class RandIter>
	inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::nonzerorandom(RandIter& g, Rep& r, Degree d) const
	{
		return random(g, r,d);
	}

	template <class Domain> template<class RandIter>
	inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::nonzerorandom(RandIter& g, Rep& r, const Rep& b) const
	{
		return random(g, r,b);
	}




	template <class Domain>
	inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::reverse( Rep& P, const Rep& Q) const {

		//     this->setdegree(Q);
		P.resize(Q.size());
		std::reverse_copy(Q.begin(), Q.end(), P.begin());
		this->setdegree(P);
		return P;
	}



	template <class Domain>
	inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::reversein( Rep& P) const
	{
		this->setdegree(P);
		std::reverse(P.begin(), P.end());
		this->setdegree(P);
		return P;
	}

} // Givaro

#endif // __GIVARO_poly_misc_INL

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
