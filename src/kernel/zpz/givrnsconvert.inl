// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givrnsconvert.inl,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givrnsconvert.inl,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description:

  // -- Computation of a mixed-radix representation of the residu.
template<class RING, class MODULO>
void RNSsystem<RING,MODULO>::RnsToMixedRadix
  (RNSsystem<RING,MODULO>::array& mixrad, const RNSsystem<RING,MODULO>::array& residu) const
{
  long i,j;
  size_t size = _primes.size();
  if (mixrad.size() < size) mixrad.reallocate( size );
  
  // -- Computation of  Ck 
  if (_ck.size()==0) ((RNSsystem*)this)->ComputeCk();

  // -- size-1 steps
  MODULO tmp;
  mixrad[0] = residu[0];
  cout << "mixrad0 = " << mixrad[0] << endl;
  for (i=1; i < size; i++)
  {  // - computes pp_i = r_0 + r_1*p_0 + ... + r_{i-1} \prod_{j<i-2} p_j [p_i]
     // Horner scheme
     tmp = mixrad[i-1]; cout << "mod" << _primes[i] << "::  " << tmp << endl;
     for (j= i-2; j>=0; j--) {
       cout << "mod" << _primes[i] << "::  " << tmp << " * " << _primes[j] << " + " << mixrad[j] <<" =";
       tmp = muladd( _primes[i], tmp, _primes[j], mixrad[j]);
       cout << " " << tmp << endl; 
     }
     cout << "\nmod" << _primes[i] << ", sumprod = " << tmp << endl;
     // - m_i = (r_i - pp_i)*ck_i, ck is reciprocals
     tmp = sub(_primes[i], residu[i], tmp);
     cout << "mod" << _primes[i] << ", sum - sumprod = " << tmp << endl;
     mixrad[i] = mul( _primes[i], tmp, _ck[i] );
     cout << "mod" << _primes[i] << ", mixrad = " << mixrad[i] << endl;
  }
}
 
  // -- Convert a mixed radix representation to an Integer
template<class RING, class MODULO>
void RNSsystem<RING,MODULO>::MixedRadixToRing( RING& res, const RNSsystem<RING,MODULO>::array& mixrad ) const 
{
  size_t size = _primes.size();
  if (size != mixrad.size()) 
    throw GivError("[RNSsystem::MixedRadixToRing]: bad size of input array");
  res = mixrad[size-1];
  for (int i=size-2; i>=0; --i) {
    res *= _primes[i];
    res += mixrad[i];
  }
}


  // Convert an integer to a RNS representation (which is given by this)
template<class RING, class MODULO>
void RNSsystem<RING,MODULO>::RingToRns( RNSsystem<RING,MODULO>::array& rns , const RING& a) const
{
  size_t size = _primes.size();
  if (rns.size() != size) rns.reallocate(size);
  // -- may be faster using the recursive 
  // tree algorithm a mod p_1...p_k/2, and a mod p_k/2+1...p_k
  for (int i=0; i<size; i++) 
    rns[i] = mod(a, _primes[i]);
}

  // Convert to an Integer:
template<class RING, class MODULO>
void RNSsystem<RING,MODULO>::RnsToRing( RING& I, const Array0<MODULO>& rns) const 
{
  // - Computation of a mixed radix representation of this
  RNSsystem<RING,MODULO>::array mixrad(_primes.size());
  RnsToMixedRadix( mixrad , rns );

  // - Convert mixrad to an integer
  MixedRadixToRing( I, mixrad ) ;
  return;
}

