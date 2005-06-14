  // -- Computation of a mixed-radix representation of the residu.
#ifndef __ECC
template<template<class> class Container> template<class TT>
inline void IntRNSsystem< Container >::RnsToMixedRadix
  (IntRNSsystem< Container >::array& mixrad, const Container<TT>& residu) const
#else
template<class Container> template<class ContTT>
inline void IntRNSsystem< Container >::RnsToMixedRadix
  (IntRNSsystem< Container >::array& mixrad, const ContTT& residu) const
#endif
{
  size_t size = _primes.size();
  if (mixrad.size() < size) mixrad.resize( size );
  
  // -- Computation of  Ck 
  if (_ck.size()==0) ((IntRNSsystem*)this)->ComputeCk();

  // -- size-1 steps
  Element tmp;
  mixrad[0] = residu[0];
  for (unsigned long i=1; i < size; i++)
  {  // - computes pp_i = r_0 + r_1*p_0 + ... + r_{i-1} \prod_{j<i-2} p_j [p_i]
     // Horner scheme
     tmp = mixrad[i-1];
// JGD 16.04.2003, Si i==1 !!!!!!!
     for (long j= i-2; j>=0; --j) {
//  std::cerr << tmp << " * " << _primes[j] << " + " << mixrad[j] << " mod " << _primes[i] << " = ";
         modin( addin( mulin(tmp, _primes[j]), mixrad[j]), _primes[i]);
//  std::cerr << tmp << ";#Horner scheme" << std::endl;
     }
     // - m_i = (r_i - pp_i)*ck_i, ck is reciprocals
//  std::cerr << "(" << residu[i] << " - " << tmp << ") * " << _ck[i] << " mod " << _primes[i] << " = ";
     mod(mixrad[i],mulin(sub(tmp,  residu[i], tmp),_ck[i]) , _primes[i] );
//  std::cerr << mixrad[i] << ";#mixrad" << std::endl;
  }
  
}
 
  

  // -- Convert a mixed radix representation to an Integer
#ifndef __ECC
template<template<class> class Container>
#else
template<class Container>
#endif
inline void IntRNSsystem< Container >::MixedRadixToRing( Element& res, const IntRNSsystem< Container >::array& mixrad ) const 
{
  size_t size = _primes.size();
//  if (size != mixrad.size()) throw GivError("[IntRNSsystem::MixedRadixToRing]: bad size of input array");
  res = mixrad[size-1];
  for (int i=size-2; i>=0; --i) {
      addin( mulin(res, _primes[i]), mixrad[i]);
//     res *= _primes[i];
//     res += mixrad[i];
  }
}


  // Convert an integer to a RNS representation (which is given by this)
#ifndef __ECC
template<template<class> class Container>
#else
template<class Container>
#endif
inline void IntRNSsystem< Container >::RingToRns( IntRNSsystem< Container >::array& rns , const external& a) const
{
  size_t size = _primes.size();
  if (rns.size() != size) rns.resize(size);
  // -- may be faster using the recursive 
  // tree algorithm a mod p_1...p_k/2, and a mod p_k/2+1...p_k
  for (int i=0; i<size; i++) 
      mod( rns[i], a, _primes[i]);
//     rns[i] = mod(a, _primes[i]);
}

  // Convert to an Integer:
#ifndef __ECC
template<template<class> class Container> template<class TT>
inline void IntRNSsystem< Container >::RnsToRing( external& I, const Container<TT>& rns) const 
#else
template<class Container> template<class ContTT>
inline void IntRNSsystem< Container >::RnsToRing( external& I, const ContTT& rns) const 
#endif
{
  // - Computation of a mixed radix representation of this
    
  typename IntRNSsystem< Container >::array mixrad(_primes.size());
  RnsToMixedRadix( mixrad , rns );

  // - Convert mixrad to an integer
  MixedRadixToRing( I, mixrad ) ;
  return;
}

