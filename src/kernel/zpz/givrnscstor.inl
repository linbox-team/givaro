// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givrnscstor.inl,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givrnscstor.inl,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description:

// -- free memory allocated in array !
template<class RING, class MODULO>
RNSsystem<RING,MODULO>::~RNSsystem()
{}

template<class RING, class MODULO>
RNSsystem<RING,MODULO>::RNSsystem ()
 : _primes(0), _ck(0)
{}

template<class RING, class MODULO>
RNSsystem<RING,MODULO>::RNSsystem (const RNSsystem<RING,MODULO>& R)
 : _primes(R._primes, givWithCopy()), 
   _ck(R._primes, givWithCopy()), 
{}


  // -- Array of primes are given
template<class RING, class MODULO>
RNSsystem<RING,MODULO>::RNSsystem( const RNSsystem<RING,MODULO>::array& inprimes) 
 : _primes(inprimes, givWithCopy()),
   _ck(0)
{
   GIVARO_ASSERT( inprimes.size()>0, "[RNSsystem<RING,MODULO>::RNSsystem] bad size of array");
}

  // -- Computes Ck , Ck = (\prod_{i=0}^{k-1} primes[i])^(-1) % primes[k],
  // for k=1..(_sz-1)
template<class RING, class MODULO>
void RNSsystem<RING,MODULO>::ComputeCk()
{
  if (_ck.size() !=0) return; // -- already computed

  // - reallocation of a new array :
  size_t size = _primes.size();
  _ck.reallocate(size);
  _ck[0] = Neutral::zero; // -- undefined and never used

  for (size_t k=1; k < size; k++)
  {
    MODULO prod = _primes[0];
    for (MODULO i= 1; i < k; i++)
       prod = mul(_primes[k], prod, _primes[i]);
    _ck[k] = inv(_primes[k], prod);
  }
}


template<class RING, class MODULO>
const RNSsystem<RING,MODULO>::array& RNSsystem<RING,MODULO>::Primes() const
{ 
  return _primes; 
}


template<class RING, class MODULO>
const MODULO RNSsystem<RING,MODULO>::ith(const size_t i) const
{
  return _primes[i];
}


template<class RING, class MODULO>
const RNSsystem<RING,MODULO>::array& RNSsystem<RING,MODULO>::Reciprocals() const
{
  if (_ck.size() ==0) ((RNSsystem<RING,MODULO>*)this)->ComputeCk();
  return _ck;
}


template<class RING, class MODULO>
const MODULO RNSsystem<RING,MODULO>::reciprocal(const size_t i) const
{
  if (_ck.size() ==0) ((RNSsystem<RING,MODULO>*)this)->ComputeCk();
  return _ck[i];
}
