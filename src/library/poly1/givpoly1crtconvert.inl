// ==========================================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: J-G Dumas
// Time-stamp: <20 May 10 15:52:06 Jean-Guillaume.Dumas@imag.fr>
// Description: Polynomial Chinese Remaindering of degree 1
// ==========================================================================
#ifndef __GIVARO_poly1_crt_convert_INL
#define __GIVARO_poly1_crt_convert_INL

namespace Givaro {
    // Converts a Polynomial to a RNS representation (which is given by this)
    template<class Field>
    typename Poly1CRT<Field>::array_T& Poly1CRT<Field>::RingToRns( typename Poly1CRT<Field>::array_T& rns , const typename Poly1CRT<Field>::Element& a) const {
        size_t Size = _primes.size();
        if (rns.size() != Size) rns.resize(Size);
        for (size_t i=0; i<Size; i++)
            _PolRing.eval(rns[i], a, _primes[i]);
        return rns;
    }
} // Givaro

#ifdef  GIVARO_CRT_EARLY_TERMINATION
#ifndef GIVARO_CRT_EARLY_TERMINATION_THRESHOLD
#define GIVARO_CRT_EARLY_TERMINATION_THRESHOLD 4
#endif
#endif


namespace Givaro {

    // Converts an array of field residues to a Polynomial
    template<class Field>
    typename Poly1CRT<Field>::Element& Poly1CRT<Field>::RnsToRing(typename Poly1CRT<Field>::Element& I, const typename Poly1CRT<Field>::array_T& rns) {
        this->ComputeCk();
        size_t Size = _primes.size();
        _PolRing.assign(I, Degree(0), rns[0]);
#ifdef  GIVARO_CRT_EARLY_TERMINATION
        size_t EarlyTerm = 0;
#endif
        // _PolRing.write(std::cerr<< "R[0]: ", I) <<std::endl;
        for(size_t i=1; i<Size; ++i) {
            Type_t addon;
            _PolRing.eval(addon, I, _primes[i]);

            //         _F.write(_PolRing.write(_F.write(std::cout << "subs(X=", _primes[i]) << ",", I) << ") mod " << _F.characteristic() << " = ", addon) << ';' << std::endl;

            _F.negin(addon);
            _F.addin(addon, rns[i]);
#ifdef  GIVARO_CRT_EARLY_TERMINATION
            if (_F.isZero(addon))
                ++EarlyTerm;
            else {
                EarlyTerm=0;
#endif // GIVARO_CRT_EARLY_TERMINATION
                _PolRing.axpyin(I, addon, _ck[i]);
                // _PolRing.modin(I, _ck[i+1]);
#ifdef  GIVARO_CRT_EARLY_TERMINATION
            }
            if (EarlyTerm >= EarlyTermThreshold)
                break;
#endif // GIVARO_CRT_EARLY_TERMINATION
        }
        //     loc.stop(); std::cerr << "RnsToRing. CRT: " << loc << std::endl;
        return I;
    }

} // Givaro
#endif // __GIVARO_poly1_crt_convert_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
