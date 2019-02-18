// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

/*!@internal
 * @file integer/givintrns_cstor.inl
 * @brief NO DOC
 */

#ifndef __GIVARO_intrns_cstor_INL
#define __GIVARO_intrns_cstor_INL

namespace Givaro {

    // -- Array of primes are given
    //#ifndef __ECC
    //template<template<class> class Container>
    //#else
    template<template <class,class> class Container, template <class> class Alloc>
    //#endif
    inline IntRNSsystem<Container, Alloc>::IntRNSsystem( const array& inprimes) :
        _primes(inprimes),
        _prod(one), _ck(0)
    {
        GIVARO_ASSERT( inprimes.size()>0, "[IntRNSsystem::IntRNSsystem] bad size of array");
    }

    // -- Array of primes are given
    //#ifndef __ECC
    template<template<class,class> class Container, template <class> class Alloc>
    template <class TT>
    inline IntRNSsystem< Container, Alloc >::IntRNSsystem( const Container<TT, Alloc<TT> >& inprimes) :
        _prod(one), _ck(0)
    {
        GIVARO_ASSERT( inprimes.size()>0, "[IntRNSsystem::IntRNSsystem] bad size of array");
        _primes.resize(inprimes.size());
        typename Container<TT, Alloc<TT> >::const_iterator np = inprimes.begin();
        for(typename array::iterator pi = _primes.begin(); pi != _primes.end(); ++pi, ++np)
            *pi = Element( *np );
    }

#if 0
    // #else
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
    // #endif
#endif

    // -- Product of primes
    //#ifndef __ECC
    template<template<class,class> class Container, template <class> class Alloc>
    //#else
    //template<class Container>
    //#endif
    inline void IntRNSsystem< Container, Alloc >::ComputeProd()
    {
        if (isOne(_prod))
            for (typename array::const_iterator pi = _primes.begin();pi != _primes.end(); ++pi)
                mulin(_prod, *pi);
    }

    // -- Computes Ck , Ck = (\prod_{i=0}^{k-1} primes[i])^(-1) % primes[k],
    // for k=1..(_sz-1)
    //#ifndef __ECC
    template<template<class,class> class Container, template <class> class Alloc>
    //#else
    //template<class Container>
    //#endif
    inline void IntRNSsystem< Container, Alloc >::ComputeCk()
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


    //#ifndef __ECC
    template<template<class,class> class Container, template<class> class Alloc>
    //#else
    //template<class Container>
    //#endif
    inline const typename IntRNSsystem< Container, Alloc >::array& IntRNSsystem< Container, Alloc >::Primes() const
    {
        return _primes;
    }


    //#ifndef __ECC
    template<template<class, class> class Container, template <class> class Alloc>
    //#else
    //template<class Container>
    //#endif
    inline const typename IntRNSsystem< Container, Alloc >::Element IntRNSsystem< Container, Alloc >::ith(const size_t i) const
    {
        return _primes[i];
    }


    //#ifndef __ECC
    template<template<class,class> class Container, template <class> class Alloc>
    //#else
    //template<class Container>
    //#endif
    inline const typename IntRNSsystem< Container, Alloc >::Element IntRNSsystem< Container, Alloc >::product() const
    {
        ((IntRNSsystem< Container, Alloc >*)this)->ComputeProd();
        return _prod;
    }

    //#ifndef __ECC
    template<template<class,class> class Container, template <class> class Alloc>
    //#else
    //template<class Container>
    //#endif
    inline const typename IntRNSsystem< Container, Alloc >::array& IntRNSsystem< Container, Alloc >::Reciprocals() const
    {
        if (_ck.size() ==0) ((IntRNSsystem< Container, Alloc >*)this)->ComputeCk();
        return _ck;
    }


    //#ifndef __ECC
    template<template<class,class> class Container, template <class> class Alloc>
    //#else
    //template<class Container>
    //#endif
    inline const typename IntRNSsystem< Container, Alloc >::Element IntRNSsystem< Container, Alloc >::reciprocal(const size_t i) const
    {
        if (_ck.size() ==0) ((IntRNSsystem< Container, Alloc >*)this)->ComputeCk();
        return _ck[i];
    }

} // namespace Givaro


#endif // __GIVARO_intrns_cstor_INL

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
