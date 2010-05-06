// ==========================================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software. 
// see the COPYRIGHT file for more details.
// Authors: J-G Dumas
// Time-stamp: <06 May 10 13:51:41 Jean-Guillaume.Dumas@imag.fr> 
// Description: Polynomial Chinese Remaindering of degree 1
// ==========================================================================

  // Converts a Polynomial to a RNS representation (which is given by this)
template<class Field>
typename Poly1CRT<Field>::array_T& Poly1CRT<Field>::RingToRns( typename Poly1CRT<Field>::array_T& rns , const typename Poly1CRT<Field>::Element& a) const {
    size_t size = _primes.size();
    if (rns.size() != size) rns.resize(size);
    for (size_t i=0; i<size; i++) 
        _PolRing.eval(rns[i], a, _primes[i]);
    return rns;
}

#ifdef  GIVARO_CRT_EARLY_TERMINATION
#ifndef GIVARO_CRT_EARLY_TERMINATION_THRESHOLD
#define GIVARO_CRT_EARLY_TERMINATION_THRESHOLD 4
#endif
#endif

  // Converts an array of field residues to a Polynomial
template<class Field>
typename Poly1CRT<Field>::Element& Poly1CRT<Field>::RnsToRing(typename Poly1CRT<Field>::Element& I, const typename Poly1CRT<Field>::array_T& rns) {
    this->ComputeCk();
    size_t size = _primes.size();
    _PolRing.assign(I, Degree(0), rns[0]);
#ifdef  GIVARO_CRT_EARLY_TERMINATION
    size_t EarlyTerm = 0;
#endif
//     _PolRing.write(std::cerr<< "R[0]: ", I) <<std::endl;
    for(size_t i=1; i<size; ++i) {
        Type_t addon; 
        _PolRing.eval(addon, I, _primes[i]);
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

