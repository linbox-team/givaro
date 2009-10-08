// ==========================================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software. 
// see the COPYRIGHT file for more details.
// Authors: J-G Dumas
// Time-stamp: <08 Oct 09 12:37:23 Jean-Guillaume.Dumas@imag.fr> 
// Description: Polynomial Chinese Remaindering of degree 1
// ==========================================================================

  // Convert an integer to a RNS representation (which is given by this)
template<class Field>
typename Poly1CRT<Field>::array_T& Poly1CRT<Field>::RingToRns( typename Poly1CRT<Field>::array_T& rns , const typename Poly1CRT<Field>::Element& a) const {
    size_t size = _primes.size();
    if (rns.size() != size) rns.reallocate(size);
    for (size_t i=0; i<size; i++) 
        _PolRing.eval(rns[i], a, _primes[i]);
    return rns;
}
#define EarlyTermThreshold 4

  // Convert to an Integer:
template<class Field>
typename Poly1CRT<Field>::Element& Poly1CRT<Field>::RnsToRing(typename Poly1CRT<Field>::Element& I, const typename Poly1CRT<Field>::array_T& rns) {
//     Timer loc; loc.start();
    this->ComputeCk();
//     loc.stop(); std::cerr << "RnsToRing. CompCK: " << loc << std::endl;
//     loc.start();
    size_t size = _primes.size();
    _PolRing.assign(I, Degree(0), rns[0]);
    size_t EarlyTerm = 0;
//     _PolRing.write(std::cerr<< "R[0]: ", I) <<std::endl;
    for(size_t i=1; i<size; ++i) {
        Type_t addon; 
        _PolRing.eval(addon, I, _primes[i]);
        _F.negin(addon);
        _F.addin(addon, rns[i]);
        if (_F.isZero(addon))
            ++EarlyTerm;
        else {
            EarlyTerm=0;
            _PolRing.axpyin(I, addon, _ck[i]);
                // _PolRing.modin(I, _ck[i+1]);
        }
        if (EarlyTerm >= EarlyTermThreshold)
            break;
    }
//     loc.stop(); std::cerr << "RnsToRing. CRT: " << loc << std::endl;
    return I;
}

