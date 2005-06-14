  // -- Array of primes are given
#ifndef __ECC
template<template<class> class Container>
#else
template<class Container>
#endif
inline IntRNSsystem< Container >::IntRNSsystem( const IntRNSsystem< Container >::array& inprimes) 
 : _primes(inprimes),
   _prod(one), _ck(0)
{
   GIVARO_ASSERT( inprimes.size()>0, "[IntRNSsystem< Container >::IntRNSsystem] bad size of array");
}

  // -- Array of primes are given
#ifndef __ECC
template<template<class> class Container> template <class TT>
inline IntRNSsystem< Container >::IntRNSsystem( const Container<TT>& inprimes) 
 : _prod(one), _ck(0)
{
   GIVARO_ASSERT( inprimes.size()>0, "[IntRNSsystem< Container >::IntRNSsystem] bad size of array");
   _primes.resize(inprimes.size());
   typename Container<TT>::const_iterator np = inprimes.begin();
   for(typename array::iterator pi = _primes.begin(); pi != _primes.end(); ++pi, ++np)
       *pi = Element( *np );
}
#else
template<class Container> template <class ContTT>
inline IntRNSsystem< Container >::IntRNSsystem( const ContTT& inprimes) 
 : _prod(one), _ck(0)
{
   GIVARO_ASSERT( inprimes.size()>0, "[IntRNSsystem< Container >::IntRNSsystem] bad size of array");
   _primes.resize(inprimes.size());
   typename ContTT::const_iterator np = inprimes.begin();
   for(typename array::iterator pi = _primes.begin(); pi != _primes.end(); ++pi, ++np)
       *pi = Element( *np );
}
#endif

  // -- Product of primes
#ifndef __ECC
template<template<class> class Container>
#else
template<class Container>
#endif
inline void IntRNSsystem< Container >::ComputeProd()
{
    if (isOne(_prod))
        for (typename array::const_iterator pi = _primes.begin();pi != _primes.end(); ++pi)
            mulin(_prod, *pi);
}

  // -- Computes Ck , Ck = (\prod_{i=0}^{k-1} primes[i])^(-1) % primes[k],
  // for k=1..(_sz-1)
#ifndef __ECC
template<template<class> class Container>
#else
template<class Container>
#endif
inline void IntRNSsystem< Container >::ComputeCk()
{
  if (_ck.size() !=0) return; // -- already computed

  // - reallocation of a new array :
  size_t size = _primes.size();
  _ck.resize(size);
  _ck[0] = Element::zero; // -- undefined and never used

  for (size_t k=1; k < size; ++k)
  {
    Element prod = _primes[0];
    for (size_t i= 1; i < k; ++i)
        modin( mulin( prod, _primes[i]), _primes[k]);
    
    Element g,u;
    gcd(g,u, _ck[k],_primes[k],prod); // _ck[k] * prod = g mod _primes[k]
  }
}


#ifndef __ECC
template<template<class> class Container>
#else
template<class Container>
#endif
inline const typename IntRNSsystem< Container >::array& IntRNSsystem< Container >::Primes() const
{ 
  return _primes; 
}


#ifndef __ECC
template<template<class> class Container>
#else
template<class Container>
#endif
inline const typename IntRNSsystem< Container >::Element IntRNSsystem< Container >::ith(const size_t i) const
{
  return _primes[i];
}


#ifndef __ECC
template<template<class> class Container>
#else
template<class Container>
#endif
inline const typename IntRNSsystem< Container >::Element IntRNSsystem< Container >::product() const
{
    ((IntRNSsystem< Container >*)this)->ComputeProd();
    return _prod;
}

#ifndef __ECC
template<template<class> class Container>
#else
template<class Container>
#endif
inline const typename IntRNSsystem< Container >::array& IntRNSsystem< Container >::Reciprocals() const
{
  if (_ck.size() ==0) ((IntRNSsystem< Container >*)this)->ComputeCk();
  return _ck;
}

 
#ifndef __ECC
template<template<class> class Container>
#else
template<class Container>
#endif
inline const typename IntRNSsystem< Container >::Element IntRNSsystem< Container >::reciprocal(const size_t i) const
{
  if (_ck.size() ==0) ((IntRNSsystem< Container >*)this)->ComputeCk();
  return _ck[i];
}
