// ==========================================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: J-G Dumas
// Time-stamp: <09 Jun 10 18:15:58 Jean-Guillaume.Dumas@imag.fr>
// Description: Polynomial Chinese Remaindering of degree 1
// ==========================================================================
#ifndef __GIVARO_poly1_crt_cstor_INL
#define __GIVARO_poly1_crt_cstor_INL

namespace Givaro {
    // -- free memory allocated in array !
    template<class Field>
    Poly1CRT<Field>::~Poly1CRT()
    {}

    template<class Field>
    Poly1CRT<Field>::Poly1CRT ()
    : _XIndet(), _F(), _PolRing(), _primes(0), _ck(0)
    {}

    template<class Field>
    Poly1CRT<Field>::Poly1CRT (const Self_t& R)
    : _XIndet(R._XIndet),
    _F(R._F),
    _PolRing(R._PolRing),
    _primes(R._primes),
    _ck(R._ck)
    {}


    // -- Array of primes are given
    template<class Field>
    Poly1CRT<Field>::Poly1CRT( const Field& F, const array_T& inprimes, const Indeter& X)
    : _XIndet(X),
    _F(F),
    _PolRing(F,X),
    _primes(inprimes),
    _ck(0)
    {
        GIVARO_ASSERT( inprimes.size()>0, "[Poly1CRT<Field>::Poly1CRT] bad size of array");
        //    for(typename array_T::const_iterator it=_primes.begin(); it!=_primes.end();++it)
        //        _F.write(std::cout, *it) << std::endl;

    }

    // -- Computes Ck , Ck = (\prod_{i=0}^{k-1} primes[i])^(-1) % primes[k],
    // for k=1..(_sz-1)
    template<class Field>
    void Poly1CRT<Field>::ComputeCk()
    {
        if (_ck.size() !=0) return; // -- already computed

        size_t Size = _primes.size();
        _ck.resize(Size+1);
        Element irred; _PolRing.init(irred, Degree(1));
        Element prod; _PolRing.init(prod, Degree(0));
        // Never used
        //   _PolRing.assign(_ck[0], prod);
        for (size_t k=1; k < Size; ++k) {
            _F.assign(irred[0], _primes[k-1]);
            _F.negin(irred[0]);
            //       _PolRing.write(std::cerr<< "irred["<<k<<"]: ", irred) <<std::endl;
            _PolRing.mulin(prod, irred);
            //       _PolRing.write(std::cerr<< "prod["<<k<<"]: ", prod) <<std::endl;
            Type_t invC; _F.init(invC);
            _PolRing.eval(invC, prod, _primes[k]);
            //       _F.write(std::cerr<< "eval["<<k<<"]: ", invC) <<std::endl;
            _F.invin(invC);
            //       _F.write(std::cerr<< "inv["<<k<<"]: ", invC) <<std::endl;
            _PolRing.mul(_ck[k],prod,invC);
            //       _PolRing.write(std::cerr<< "mul["<<k<<"]: ", _ck[k]) <<std::endl;
        }
        _F.assign(irred[0], _primes[Size-1]);
        _F.negin(irred[0]);
        _PolRing.mul(_ck[Size], prod, irred);

        //    for(typename array_E::const_iterator it=_ck.begin(); it!=_ck.end();++it)
        //        _PolRing.write(std::cout, *it) << std::endl;

    }


    template<class Field>
    const typename Poly1CRT<Field>::array_T& Poly1CRT<Field>::Primes() const
    {
        return _primes;
    }


    template<class Field>
    const typename Poly1CRT<Field>::Type_t& Poly1CRT<Field>::ith(const size_t i) const
    {
        return _primes[i];
    }


    template<class Field>
    const typename Poly1CRT<Field>::array_E& Poly1CRT<Field>::Reciprocals() const
    {
        if (_ck.size() ==0) ((Poly1CRT<Field>*)this)->ComputeCk();
        return _ck;
    }


    template<class Field>
    const typename Poly1CRT<Field>::Element& Poly1CRT<Field>::reciprocal(const size_t i) const
    {
        if (_ck.size() ==0) ((Poly1CRT<Field>*)this)->ComputeCk();
        return _ck[i];
    }
} // Givaro
#endif // __GIVARO_poly1_crt_cstor_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
