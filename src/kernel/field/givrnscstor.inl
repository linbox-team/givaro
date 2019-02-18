// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givrnscstor.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gaut%ier
// $Id: givrnscstor.inl,v 1.4 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================
// Description:

#ifndef __GIVARO_rns_cstor_INL
#define __GIVARO_rns_cstor_INL

namespace Givaro {

	// -- free memory allocated in array !
	template<class RING, class Domain>
	RNSsystem<RING,Domain>::~RNSsystem()
	{
	}

	template<class RING, class Domain>
	RNSsystem<RING,Domain>::RNSsystem ()
	: _primes(0), _ck(0)
	{}

	template<class RING, class Domain>
	RNSsystem<RING,Domain>::RNSsystem (const Self_t & R)
	: _primes(R._primes, givWithCopy()),
	_ck(R._ck, givWithCopy())
	{}

	template<class RING, class Domain>
	void RNSsystem<RING,Domain>::setPrimes (const domains& inprimes)
	{
		_primes.allocate(0);
		_primes.copy( inprimes );
		_ck.resize(0);
	}

	// -- Array of primes are given
	template<class RING, class Domain>
	RNSsystem<RING,Domain>::RNSsystem( const domains& inprimes)
	: _primes(inprimes, givWithCopy()),
	_ck(0)
	{
		GIVARO_ASSERT( inprimes.size()>0, "[RNSsystem<RING,Domain>::RNSsystem] bad size of array");
	}

	// -- Computes Ck , Ck = (\prod_{i=0}^{k-1} primes[i])^(-1) % primes[k],
	// for k=1..(_sz-1)
	template<class RING, class Domain>
	void RNSsystem<RING,Domain>::ComputeCk()
	{
		if (_ck.size() !=0) return; // -- already computed

		// - reallocation of a new array :
		int Size = (int) _primes.size();
		_ck.resize((size_t)Size);
		//  _ck[0] = Neutral::zero; // -- undefined and never used

		for (int k=1; k < Size; ++k)
		{
			modulo prod, tmp;
			_primes[k].init(prod, _primes[0].characteristic());
			for (int i= 1; i < k; ++i)
				_primes[k].mulin(prod, _primes[k].init(tmp,_primes[i].characteristic()));
			_primes[k].inv(_ck[k],  prod);
		}
	}


	template<class RING, class Domain>
	const typename RNSsystem<RING,Domain>::domains& RNSsystem<RING,Domain>::Primes() const
	{
		return _primes;
	}


	template<class RING, class Domain>
	const Domain RNSsystem<RING,Domain>::ith(const size_t i) const
	{
		return _primes[i];
	}


	template<class RING, class Domain>
	const typename RNSsystem<RING,Domain>::array& RNSsystem<RING,Domain>::Reciprocals() const
	{
		if (_ck.size() ==0) ((RNSsystem<RING,Domain>*)this)->ComputeCk();
		return _ck;
	}


	template<class RING, class Domain>
	const typename RNSsystem<RING,Domain>::modulo RNSsystem<RING,Domain>::reciprocal(const size_t i) const
	{
		if (_ck.size() ==0) ((RNSsystem<RING,Domain>*)this)->ComputeCk();
		return _ck[i];
	}

} // namespace Givaro

#endif // __GIVARO_rns_cstor_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
