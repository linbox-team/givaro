// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

#ifndef __GIVARO_intrns_convert_INL
#define __GIVARO_intrns_convert_INL

namespace Givaro {

    // -- Computation of a mixed-radix representation of the residu.
    //#ifndef __ECC
    template<template<class, class> class Container, template <class> class Alloc>
    template<class TT>
    inline void IntRNSsystem< Container, Alloc >::RnsToMixedRadix
    (IntRNSsystem< Container, Alloc >::array& mixrad, const Container<TT, Alloc<TT> >& residu)
    //#else
    //template<class Container> template<class ContTT>
    //inline void IntRNSsystem< Container >::RnsToMixedRadix
    //  (IntRNSsystem< Container >::array& mixrad, const ContTT& residu) const
    //#endif
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
            for (long j= (long)i-1; j--; ) {
                //  std::cerr << tmp << " * " << _primes[j] << " + " << mixrad[j] << " mod " << _primes[i] << " = ";
                modin( addin( mulin(tmp, _primes[(size_t)j]), mixrad[(size_t)j]), _primes[i]);
                //  std::cerr << tmp << ";#Horner scheme" << std::endl;
            }
            // - m_i = (r_i - pp_i)*ck_i, ck is reciprocals
            //  std::cerr << "(" << residu[i] << " - " << tmp << ") * " << _ck[i] << " mod " << _primes[i] << " = ";
            mod(mixrad[i],mulin(sub(tmp,  residu[i], tmp),_ck[i]) , _primes[i] );
            //  std::cerr << mixrad[i] << ";#mixrad" << std::endl;
        }

    }



    // -- Convert a mixed radix representation to an Integer
    //#ifndef __ECC
    template<template<class,class> class Container, template <class> class Alloc>
    //#else
    //template<class Container>
    //#endif
    inline void IntRNSsystem< Container, Alloc >::MixedRadixToRing( Element& res, const IntRNSsystem< Container, Alloc >::array& mixrad ) const
    {
        size_t size = _primes.size();
        //  if (size != mixrad.size()) throw GivError("[IntRNSsystem::MixedRadixToRing]: bad size of input array");
        res = mixrad[size-1];
        for (size_t i=size-1; i--; ) {
            addin( mulin(res, _primes[i]), mixrad[i]);
            //     res *= _primes[i];
            //     res += mixrad[i];
        }
    }


    // Convert an integer to a RNS representation (which is given by this)
    //#ifndef __ECC
    template<template<class,class> class Container,template<class> class Alloc>
    //#else
    //template<class Container>
    //#endif
    inline void IntRNSsystem< Container, Alloc >::RingToRns( IntRNSsystem< Container, Alloc >::array& rns , const external& a)
    {
        size_t mysize = _primes.size();
        if (rns.size() != mysize) rns.resize(mysize);
        // -- may be faster using the recursive
        // tree algorithm a mod p_1...p_k/2, and a mod p_k/2+1...p_k
        for (int i=0; i<(int)mysize; i++)
            mod( rns[i], a, _primes[i]);
        //     rns[i] = mod(a, _primes[i]);
    }

    // Convert to an Integer:
    //#ifndef __ECC
    template<template<class, class> class Container, template <class> class Alloc>
    template<class TT>
    inline void IntRNSsystem< Container, Alloc >::RnsToRing( external& I, const Container<TT, Alloc<TT> >& rns)
    //#else
    //template<class Container> template<class ContTT>
    //inline void IntRNSsystem< Container >::RnsToRing( external& I, const ContTT& rns) const
    //#endif
    {
        // - Computation of a mixed radix representation of this

        typename IntRNSsystem< Container, Alloc >::array mixrad(_primes.size());
        RnsToMixedRadix( mixrad , rns );

        // - Convert mixrad to an integer
        MixedRadixToRing( I, mixrad ) ;
        return;
    }

} // namespace Givaro

#endif // __GIVARO_intrns_convert_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
