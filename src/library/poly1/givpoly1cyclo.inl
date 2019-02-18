// =============================================================== //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Givaro / Athapascan-1
// Cyclotomic polynomials
// Time-stamp: <15 Jul 08 10:40:47 Jean-Guillaume.Dumas@imag.fr>
// =============================================================== //
#ifndef __GIVARO_poly1_cyclo_INL
#define __GIVARO_poly1_cyclo_INL
#include <givaro/givintfactor.h>

namespace Givaro {

    // ---------------------------------------------------------------
    // Composition by a power of X
    // ---------------------------------------------------------------
    template<class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::power_compose(Rep& W, const Rep& P, uint64_t b) const
    {
        Degree dp; degree(dp, P);
        Type_t lc;
        leadcoef(lc, P);
        assign( W, b*dp.value(), lc); // all coeffs to zero, except leading ...
        for(size_t i=0;i<dp.value();++i) {
            _domain.assign(W[i*b], P[i]);
        }
        return setdegree(W);
    }




    // ---------------------------------------------------------------
    // n th Cyclotomic polynomial
    // ---------------------------------------------------------------
    template<class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::cyclotomic( Rep& P, uint64_t n) const
    {
        // P must provide an indeterminate
        //   Integer In(n);
        IntFactorDom<> IF;
        if (n <= 1) {
            init(P, Degree(1));
            _domain.assign(P[0], _domain.one);
            _domain.negin(P[0]);
            return setdegree(P);
        } else if (IF.isprime(n)) {
            init(P, Degree(n-1));
            for(size_t i=(size_t)n-1;--i;)
                _domain.assign(P[i], _domain.one);
            return setdegree(P);
        }
        else {
            uint64_t q,f;
            q = n / 2;
            f = n % 2;
            if (f) {
                IntFactorDom<>::Rep If, Iq;
                IF.divexact(Iq, n, IF.factor(If,n) );
                IF.convert(f, If);
                IF.convert(q, Iq);
                Rep inter;
                cyclotomic(inter,q);
                if (q % f) {
                    Rep res;
                    power_compose(res, inter, f);
                    return div(P, res, inter);
                } else
                    return power_compose(P,inter,f);
            } else {
                if (q%2) {
                    // q odd
                    cyclotomic(P,q);
                    Degree di;
                    degree(di, P);
                    for(int i=1;i<=di.value();i+=2)
                        _domain.negin(P[i]);
                    return setdegree(P);
                } else {
                    // q even
                    Rep inter;
                    cyclotomic(inter,q);
                    Degree di;
                    degree(di, inter);
                    return power_compose(P, inter, 2);
                }
            }
        }
    }

} // Givaro

#endif // __GIVARO_poly1_cyclo_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
