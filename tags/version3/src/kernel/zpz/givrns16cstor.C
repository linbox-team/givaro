// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givrns16cstor.C,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givrns16cstor.C,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description:

#include "givaro/givrns16.h"

// -- free memory allocated in array !
RNS16system::~RNS16system()
{}

RNS16system::RNS16system ()
 : _primes(0), _ck(0), _sizek(0), _log2k(0), _u(0), _qk(0)
{}

RNS16system::RNS16system (const RNS16system& R)
 : _primes(R._primes, givWithCopy()), 
   _ck(R._primes, givWithCopy()), 
   _sizek(R._sizek),
   _log2k(R._log2k),
   _u( R._sizek + R._sizek * R._log2k),
   _qk(R._qk, givWithCopy())
{}


  // -- Array of primes are given
RNS16system::RNS16system( const RNS16system::array& inprimes) 
 : _primes(inprimes, givWithCopy()),
   _ck(0), _sizek(0), _log2k(0), _u(0), _qk(0)
{
  GIVARO_ASSERT( inprimes.size()>0, "[RNS16system::RNS16system] bad size of array");
}

  // -- Computes Ck , Ck = (\prod_{i=0}^{k-1} primes[i])^(-1) % primes[k],
  // for k=1..(_sz-1)
void RNS16system::ComputeCk()
{
  if (_ck.size() !=0) return; // -- already computed

  // - reallocation of a new array :
  size_t size = _primes.size();
  _ck.reallocate(size);
  _ck[0] = Neutral::zero; // -- undefined and never used

  for (size_t k=1; k < size; k++)
  {
    RNS16system::Modulo prod = _primes[0];
    ZpzDom<Std16> Fk (_primes[k]); 
    // JGD 27.03.03
    // for (RNS16system::Modulo i= 1; i < k; i++)
    for (size_t i= 1; i < k; i++)
       Fk.mulin(prod, _primes[i]);
    Fk.inv(_ck[k], prod);
  }
}

void RNS16system::ComputeQk()
{
  if (_qk.size() !=0) return; // -- already computed

  size_t k, size_primes = _primes.size();
  _log2k = 0;
  _sizek = 1;
  k = size_primes;
  // -- compute log2k, and sizek such that 2^(sizek-1) < k <= 2^sizek
  while ( k > _sizek ) {++_log2k; _sizek <<= 1; } 

  // - reallocation of a new array :
  _qk.reallocate( _sizek * _log2k );
  _u.reallocate( _sizek + _sizek * _log2k );

  size_t i,j;
  for (i=0; i<size_primes; i++) _qk[i] = _primes[i];
  for (i=size_primes; i<_sizek; i++) _qk[i] = 1;
  for (j=1; j<_log2k; ++j) {
    size_t step = 1<<j;
    for (i=0; i<_sizek; i += step)
      //val(_qk,i,j) = val(_qk,i,j-1) * val(_qk,i+step/2, j-1);
      Integer::mul(_qk[i+_sizek*j], _qk[i+_sizek*(j-1)], _qk[i+step/2 + _sizek*(j-1)]);
  }
}




const RNS16system::array& RNS16system::Primes() const
{ 
  return _primes; 
}


const RNS16system::Modulo RNS16system::ith(const size_t i) const
{
  return _primes[i];
}


const RNS16system::array& RNS16system::Reciprocals() const
{
  if (_ck.size() ==0) ((RNS16system*)this)->ComputeCk();
  return _ck;
}


const RNS16system::Modulo RNS16system::reciprocal(const size_t i) const
{
  if (_ck.size() ==0) ((RNS16system*)this)->ComputeCk();
  return _ck[i];
}
