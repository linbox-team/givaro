// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givrnsconvert.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givrnsconvert.inl,v 1.5 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================
// Description:

#ifndef __GIVARO_rns_convert_INL
#define __GIVARO_rns_convert_INL

namespace Givaro {

  // -- Computation of a mixed-radix representation of the residu.
template<class RING, class Domain>
void RNSsystem<RING,Domain>::RnsToMixedRadix
  (RNSsystem<RING,Domain>::array& mixrad, const RNSsystem<RING,Domain>::array& residu)
{
  int j;
  int i=0,Size = (int)_primes.size();
  if ((int)mixrad.size() < Size) mixrad.resize( (size_t)Size );

  // -- Computation of  Ck
  if (_ck.size()==0) ((RNSsystem*)this)->ComputeCk();

  // -- Size-1 steps
  modulo t2, t3, t4, tmp;
  _primes[0].assign(mixrad[0],residu[0]);
// _primes[0].write(std::cerr << "mixrad0 = ", mixrad[0]) << std::endl;
  for (i=1; i < (int)Size; ++i)
  {  // - computes pp_i = r_0 + r_1*p_0 + ... + r_{i-1} \prod_{j<i-2} p_j [p_i]
     // Horner scheme

      _primes[i].init(tmp, _primes[i-1].convert(t4, mixrad[i-1]) );
// _primes[i].write(std::cerr << "within ") << " --> ";
// _primes[i].write(std::cerr << t4 << " mod " << _primes[i].characteristic() << "::  ", tmp) << std::endl;
     for (j= i-2; j>=0; --j) {
// _primes[j].write(_primes[i].write(std::cerr << "mod " << _primes[i].characteristic() << "::  ", tmp) << " * " << _primes[j].characteristic() << " + ", mixrad[j]) <<" =";

       _primes[i].init(t3, _primes[j].convert(t4, mixrad[j]));
       _primes[i].init(t4, _primes[j].characteristic());
       _primes[i].axpy(t2, tmp, t4, t3);
       _primes[i].assign(tmp, t2);
// _primes[i].write(std::cerr << " ", tmp) << std::endl;
     }
// _primes[i].write(std::cerr << "\nmod" << _primes[i].characteristic() << ", sumprod = ", tmp) << std::endl;
     // - m_i = (r_i - pp_i)*ck_i, ck is reciprocals
     _primes[i].sub(t2, residu[i], tmp);
     _primes[i].assign(tmp, t2);
// _primes[i].write(std::cerr << "mod " << _primes[i].characteristic() << ", sum - sumprod = ", tmp) << std::endl;

     _primes[i].mul(mixrad[i], tmp, _ck[i]);
// _primes[i].write(std::cerr << "mod " << _primes[i].characteristic() << ", mixrad = ", mixrad[i]) << std::endl;
  }
}

  // -- Convert a mixed radix representation to an Integer
template<class RING, class Domain>
RING& RNSsystem<RING,Domain>::MixedRadixToRing( RING& res, const RNSsystem<RING,Domain>::array& mixrad ) const
{
  size_t Size = _primes.size();
  if (!Size)
	  GivError("_primes is empty");
  if (Size != mixrad.size())
    throw GivError("[RNSsystem::MixedRadixToRing]: bad size of input array");
  _primes[int(Size-1)].convert(res,mixrad[int(Size-1)]);
  RING tmp;
  if( Size  == 1 )
	  return res;

  for (int i=int(Size-1); i--; ) {
    res *= _primes[i].characteristic();
    res += _primes[i].convert(tmp, mixrad[i]);
  }
  return res;
}


  // Convert an integer to a RNS representation (which is given by this)
template<class RING, class Domain>
void RNSsystem<RING,Domain>::RingToRns( RNSsystem<RING,Domain>::array& rns , const RING& a) const
{
  int Size = (int) _primes.size();
  if ((int)rns.size() != Size) rns.resize((size_t)Size);
  // -- may be faster using the recursive
  // tree algorithm a mod p_1...p_k/2, and a mod p_k/2+1...p_k
  for (int i=0; i<Size; i++)
      _primes[i].init(rns[i], a);
}

  // Convert to an Integer:
template<class RING, class Domain>
RING& RNSsystem<RING,Domain>::RnsToRing( RING& I, const RNSsystem<RING,Domain>::array& rns)
{
  // - Computation of a mixed radix representation of this
  typename RNSsystem<RING,Domain>::array mixrad(_primes.size());
  RnsToMixedRadix( mixrad , rns );

  // - Convert mixrad to an integer
  return MixedRadixToRing( I, mixrad ) ;
}

} // namespace Givaro

#endif // __GIVARO_rns_convert_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
