// ==========================================================================
// Copyright(c)'1994-2022 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: B. Grenet
//
// Generalized middle product:
// For polynomials P of size m+n-1 and Q of size n, define their
// middle product MP(P,Q) as the central m coefficients of P*Q
// that is MP(P,Q) = (PQ div X^(n-1)) mod X^m
//
// Balanced case: m = n
// ==========================================================================
#ifndef __GIVARO_poly1_midmul_INL
#define __GIVARO_poly1_midmul_INL

namespace Givaro {

#ifndef KARA_THRESHOLD
#define KARA_THRESHOLD 50
#endif

    // forces standard middle product
    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::stdmidmul(
        Rep& R, const Rep& P, const Rep& Q ) const
    {
        const size_t sP = P.size();
        const size_t sQ = Q.size();
        const size_t sR = sP-sQ+1;
        if (sR != R.size()) R.resize(sR);

        stdmidmul(R, R.begin(), R.end(),
               P, P.begin(), P.end(),
               Q, Q.begin(), Q.end());

        return setdegree(R);
    }

    // Standard middle product between iterator bounds
    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::stdmidmul(
        Rep& R, const RepIterator Rbeg, const RepIterator Rend,
        const Rep& P, const RepConstIterator Pbeg, const RepConstIterator Pend,
        const Rep& Q, const RepConstIterator Qbeg, const RepConstIterator Qend )
        const {

        if (Rbeg == Rend) return R;

        RepConstIterator ai=Pbeg,aig=Pbeg,bi=Qend; --bi;
        RepIterator ri=Rbeg;
        if (_domain.isZero(*bi))
            for(; (ai != Pend) && (ri != Rend); ++ai,++ri)
                _domain.assign(*ri,_domain.zero);
        else
            for(; (ai!=Pend) && (ri!=Rend);++ai,++ri)
                if (_domain.isZero(*ai))
                    _domain.assign(*ri, _domain.zero);
                else
                    _domain.mul(*ri,*ai,*bi);

        for(;ri!=Rend;++ri)
            _domain.assign(*ri,_domain.zero);

        for(--bi,++aig; (aig!=Pend) && (bi>=Qbeg);++aig,--bi)
            if (! _domain.isZero(*bi))
                for(ri=Rbeg,ai=aig; (ai!=Pend) && (ri!=Rend);++ai,++ri)
                    _domain.axpyin(*ri,*ai,*bi);
        return R;
    }

}
#endif // __GIVARO_poly1_mid_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
