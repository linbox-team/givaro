// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givrnscstor.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software. 
// see the COPYRIGHT file for more details.
// Authors: T. Gaut%ier
// $Id: givrnscstor.inl,v 1.3 2009-09-17 14:28:23 jgdumas Exp $
// ==========================================================================
// Description:

// -- free memory allocated in array !
template<class RING, class Domain>
RNSsystem<RING,Domain>::~RNSsystem()
{}

template<class RING, class Domain>
RNSsystem<RING,Domain>::RNSsystem ()
 : _primes(0), _ck(0)
{}

template<class RING, class Domain>
RNSsystem<RING,Domain>::RNSsystem (const RNSsystem<RING,Domain>& R)
 : _primes(R._primes, givWithCopy()), 
   _ck(R._ck, givWithCopy())
{}


  // -- Array of primes are given
template<class RING, class Domain>
RNSsystem<RING,Domain>::RNSsystem( const RNSsystem<RING,Domain>::domains& inprimes) 
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
  size_t size = _primes.size();
  _ck.reallocate(size);
//  _ck[0] = Neutral::zero; // -- undefined and never used

  for (size_t k=1; k < size; ++k)
  {
    modulo prod, tmp;
    _primes[k].init(prod, _primes[0].characteristic());
    for (size_t i= 1; i < k; ++i)
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
