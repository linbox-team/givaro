#ifndef __GIVARO_POLY1_CYCLO_TABLE_H_
#define __GIVARO_POLY1_CYCLO_TABLE_H_

// ==========================================================================
// Copyright(c)'1994-2000 by Givaro Team
// see the copyright file.
// Time-stamp: <05 Apr 00 10:17:06 Jean-Guillaume.Dumas@imag.fr>
// ==========================================================================

// Tabulation of factors of cyclotomic polynomials 
// of degree expo modulo mod
// By vectors. P = v[0] + v[1] X + ... + v[n] X^n

#include "givaro/givpoly1.h"
#include "givaro/givpoly1factor.h"

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

#include "givaro/givpoly1tabcycl.inl"
    
#endif // __GIVARO_POLY1_CYCLO_TABLE_H_
