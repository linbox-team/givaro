// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givrnsconvert.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software. 
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givrnsconvert.inl,v 1.4 2009-09-17 14:28:23 jgdumas Exp $
// ==========================================================================
// Description:

  // -- Computation of a mixed-radix representation of the residu.
template<class RING, class Domain>
void RNSsystem<RING,Domain>::RnsToMixedRadix
  (RNSsystem<RING,Domain>::array& mixrad, const RNSsystem<RING,Domain>::array& residu) const
{
  long j;
  size_t i=0,size = _primes.size();
  if (mixrad.size() < size) mixrad.reallocate( size );
  
  // -- Computation of  Ck 
  if (_ck.size()==0) ((RNSsystem*)this)->ComputeCk();

  // -- size-1 steps
  modulo t2, t3, t4, tmp;
  _primes[0].assign(mixrad[0],residu[0]);
//_primes[0].write(std::cout << "mixrad0 = ", mixrad[0]) << std::endl;
  for (i=1; i < size; ++i)
  {  // - computes pp_i = r_0 + r_1*p_0 + ... + r_{i-1} \prod_{j<i-2} p_j [p_i]
     // Horner scheme
      
      _primes[i].init(tmp, _primes[i-1].convert(t4, mixrad[i-1]) );
//_primes[i].write(std::cout << "mod " << _primes[i].characteristic() << "::  ", tmp) << std::endl;
     for (j= i-2; j>=0; --j) {
//_primes[j].write(_primes[i].write(std::cout << "mod " << _primes[i].characteristic() << "::  ", tmp) << " * " << _primes[j].characteristic() << " + ", mixrad[j]) <<" =";
       
       _primes[i].init(t3, _primes[j].convert(t4, mixrad[j]));
       _primes[i].init(t4, _primes[j].characteristic());
       _primes[i].axpy(t2, tmp, t4, t3);
       _primes[i].assign(tmp, t2);
//_primes[i].write(std::cout << " ", tmp) << std::endl; 
     }
//_primes[i].write(std::cout << "\nmod" << _primes[i].characteristic() << ", sumprod = ", tmp) << std::endl;
     // - m_i = (r_i - pp_i)*ck_i, ck is reciprocals
     _primes[i].sub(t2, residu[i], tmp);
     _primes[i].assign(tmp, t2);
//_primes[i].write(std::cout << "mod " << _primes[i].characteristic() << ", sum - sumprod = ", tmp) << std::endl;
     
     _primes[i].mul(mixrad[i], tmp, _ck[i]);
//_primes[i].write(std::cout << "mod " << _primes[i].characteristic() << ", mixrad = ", mixrad[i]) << std::endl;
  }
}
 
  // -- Convert a mixed radix representation to an Integer
template<class RING, class Domain>
void RNSsystem<RING,Domain>::MixedRadixToRing( RING& res, const RNSsystem<RING,Domain>::array& mixrad ) const 
{
  size_t size = _primes.size();
  if (size != mixrad.size()) 
    throw GivError("[RNSsystem::MixedRadixToRing]: bad size of input array");
  _primes[size-1].convert(res,mixrad[size-1]);
  RING tmp;
  for (int i=size-2; i>=0; --i) {
    res *= _primes[i].characteristic();
    res += _primes[i].convert(tmp, mixrad[i]);
  }
}


  // Convert an integer to a RNS representation (which is given by this)
template<class RING, class Domain>
void RNSsystem<RING,Domain>::RingToRns( RNSsystem<RING,Domain>::array& rns , const RING& a) const
{
  size_t size = _primes.size();
  if (rns.size() != size) rns.reallocate(size);
  // -- may be faster using the recursive 
  // tree algorithm a mod p_1...p_k/2, and a mod p_k/2+1...p_k
  for (size_t i=0; i<size; i++) 
      _primes[i].init(rns[i], a);
}

  // Convert to an Integer:
template<class RING, class Domain>
void RNSsystem<RING,Domain>::RnsToRing( RING& I, const RNSsystem<RING,Domain>::array& rns) const 
{
  // - Computation of a mixed radix representation of this
  typename RNSsystem<RING,Domain>::array mixrad(_primes.size());
  RnsToMixedRadix( mixrad , rns );

  // - Convert mixrad to an integer
  MixedRadixToRing( I, mixrad ) ;
  return;
}

