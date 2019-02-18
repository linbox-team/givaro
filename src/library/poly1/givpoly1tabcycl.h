// ==========================================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <05 Apr 00 10:17:06 Jean-Guillaume.Dumas@imag.fr>
// ==========================================================================
/** @file givpoly1tabcycl.h
 * @ingroup poly1
 * @brief Tabulation of factors of cyclotomic polynomials of degree expo modulo mod
 *
 * By vectors. P = v[0] + v[1] X + ... + v[n] X^n
 */
#ifndef __GIVARO_poly1_cyclo_table_H
#define __GIVARO_poly1_cyclo_table_H

#include "givaro/givpoly1.h"
#include "givaro/givpoly1factor.h"

namespace Givaro {

    //! CyclotomicTable
    template<class Domain, class Tag>
    class CyclotomicTable : public Poly1FactorDom<Domain,Tag> {
        typedef typename Poly1FactorDom<Domain,Tag>::Element Element;
        Element _Irreductible;
    public:
        CyclotomicTable(Domain& _d, const long expo, const Indeter& X = Indeter() ) : Poly1FactorDom<Domain,Tag>(_d, X), _Irreductible(zero) {
            typename Domain::Residu_t mod = _d.residu();
            table_0(mod, expo);
        }

        Element& getcyclo(Element& res) const { return res = _Irreductible; }
        void set_random_irreducible(const Domain& _d, const long expo) {
            random_irreducible( _Irreductible, expo);
        }

        void table_0 (const typename Domain::Residu_t mod, const long expo) ;
        void table_50 (const typename Domain::Residu_t mod, const long expo) ;
    };
} // Givaro

#include "givaro/givpoly1tabcycl.inl"

#endif // __GIVARO_poly1_cyclo_table_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
