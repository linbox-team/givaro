// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givrns16convert.C,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givrns16convert.C,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description:
#include "givaro/givrns16.h"


  // -- Computation of a mixed-radix representation of the residu.
void RNS16system::RnsToMixedRadix
  (RNS16system::array& mixrad, const RNS16system::array& residu) const
{
  unsigned long i,j;
  size_t size = _primes.size();
  if (mixrad.size() < size) mixrad.reallocate( size );
  
  // -- Computation of  Ck 
  if (_ck.size()==0) ((RNS16system*)this)->ComputeCk();

  // -- size-1 steps
  RNS16system::Modulo tmp;
  mixrad[0] = residu[0];
  for (i=1; i < size; i++)
  {  
     ZpzDom<Std16> Fp(_primes[i] );
     // - computes pp_i = r_0 + r_1*p_0 + ... + r_{i-1} \prod_{j<i-2} p_j [p_i]
     // Horner scheme
     tmp = mixrad[i-1]; 
     for (j= i-2; j>=0; j--) 
       Fp.axpy( tmp, tmp, _primes[j], mixrad[j]);
     // - m_i = (r_i - pp_i)*ck_i, ck is reciprocals
     Fp.sub( tmp, residu[i], tmp);
     Fp.mul( mixrad[i], tmp, _ck[i] );
  }
}
 
  // -- Convert a mixed radix representation to an Integer
void RNS16system::MixedRadixToRing( Integer& res, const RNS16system::array& mixrad ) const 
{
  size_t size = _primes.size();
  if (size != mixrad.size()) 
    throw GivError("[RNS16system::MixedRadixToRing]: bad size of input array");
  res = mixrad[size-1];
  for (int i=size-2; i>=0; --i) {
    Integer::mulin(res, (unsigned long)_primes[i]);
    Integer::addin(res, (unsigned long)mixrad[i]);
  }
}


  // Convert an integer to a RNS16 representation (which is given by this)
void RNS16system::RingToRns( RNS16system::array& residu , const Integer& a) const
{
  size_t k_moduli = _primes.size();
  if (residu.size() != k_moduli) residu.reallocate(k_moduli);
  // -- faster algorithm: next one
  Integer res =1;
  for (unsigned int i=0; i<k_moduli; i++) {
    Integer::mod(res, a, (unsigned long)_primes[i]);
    residu[i] = res[0];
  }
}

  // Convert to an Integer:
void RNS16system::RnsToRing( Integer& I, const array& rns) const 
{
  // - Computation of a mixed radix representation of this
  RNS16system::array mixrad(_primes.size());
  RnsToMixedRadix( mixrad , rns );

  // - Convert mixrad to an integer
  MixedRadixToRing( I, mixrad ) ;
  return;
}

  // Convert an integer to a RNS16 representation (which is given by this)
  // See Aho, Hoptcroft & Ullman, The design & analysis of computer algorithms.
void RNS16system::fastRingToRns( RNS16system::array& residu , const Integer& a) const
{
  size_t k_moduli = _primes.size();
  if (residu.size() != k_moduli) residu.reallocate(k_moduli);
  if (k_moduli <=2) {
    Integer res = 1;
    for (unsigned int i=0; i<k_moduli; i++) {
      Integer::mod(res, a, (unsigned long)_primes[i]);
      residu[i] = res[0]; // -- low word of Integer
    }
    return;
  } 
  if (_log2k ==0) ((RNS16system*)this)->ComputeQk();
  size_t i,j;
  Array0<Integer>& U = (Array0<Integer>&)_u;
  U[_sizek*_log2k] = a;
  for (j=_log2k; j>=1; --j) {
    size_t step = 1<<j;
    size_t ind1 = _sizek*(j-1);
    size_t ind2 = ind1 + step/2; 
    for (i=0; i<_sizek; i +=step, ind1 +=step, ind2 +=step) {
      //val(U,i,j-1) = val(U,i,j) % val(Q,i,j-1);
      //Integer::mod(U[ind1-_sizek], U[ind1], _qk[ind1-_sizek]);
      //Integer::mod(U[i+_sizek*(j-1)], U[i+_sizek*j], _qk[i+_sizek*(j-1)]);
      Integer::mod(U[ind1], U[ind1+_sizek], _qk[ind1]);

      //val(U,i+step/2,j-1) = val(U,i,j) % val(Q,i+step/2,j-1);
      //Integer::mod(U[i+step/2+_sizek*(j-1)], U[i+_sizek*j], _qk[ind2-_sizek]);
      //Integer::mod(U[i+step/2+_sizek*(j-1)], U[i+_sizek*j], _qk[i+step/2+_sizek*(j-1)]);
      Integer::mod(U[ind2], U[ind1+_sizek], _qk[ind2]);
    }
  }
  for (i=0; i < k_moduli; ++i) residu[i] = Integer2long(U[i]);
}

