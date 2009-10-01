// ==========================================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software. 
// see the COPYRIGHT file for more details.
// Authors: J-G Dumas
// Time-stamp: <30 Sep 09 16:14:45 Jean-Guillaume.Dumas@imag.fr> 
// Description: Polynomial Chinese Remaindering of degree 1
// ==========================================================================

  // Convert an integer to a RNS representation (which is given by this)
template<class Field>
typename Poly1CRT<Field>::array_T& Poly1CRT<Field>::RingToRns( typename Poly1CRT<Field>::array_T& rns , const typename Poly1CRT<Field>::Element& a) const {
    size_t size = _primes.size();
    if (rns.size() != size) rns.reallocate(size);
    for (size_t i=0; i<size; i++) 
        _P.eval(rns[i], a, _primes[i]);
    return rns;
}

  // Convert to an Integer:
template<class Field>
typename Poly1CRT<Field>::Element& Poly1CRT<Field>::RnsToRing(typename Poly1CRT<Field>::Element& I, const typename Poly1CRT<Field>::array_T& rns) {
    this->ComputeCk();
    size_t size = _primes.size();
    _P.assign(I, Degree(0), rns[0]);
//     _P.write(std::cerr<< "R[0]: ", I) <<std::endl;
    for(size_t i=1; i<size; ++i) {
        Element Cmun; _P.sub(Cmun, _ck[i], _F.one);
//         _P.write(std::cerr<< "Cmun["<<i<<"]: ", Cmun) <<std::endl;
        Element Ct; _P.mul(Ct, _ck[i], rns[i]);
//         _P.write(std::cerr<< "Ct["<<i<<"]: ", Ct) <<std::endl;
        _P.maxpy(I, Cmun, I, Ct);
//         _P.write(std::cerr<< "R["<<i<<"]: ", I) <<std::endl;
        _P.modin(I, _ck[i+1]);
//         _P.write(std::cerr<< "R["<<i<<"]: ", I) <<std::endl;
    }
    return I;
}

