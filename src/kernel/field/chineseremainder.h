// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/chineseremainder.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: chineseremainder.h,v 1.12 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================

/*!@file chineseremainder.h
 * @ingroup zpz
 * @brief  Chinese Remainder Algorithm for 2 Elements.
 * @sa
 * For any number of moduli see zpz/givrns.h or zpz/givrnsfixed.h
 */

#ifndef __GIVARO_cra_H
#define __GIVARO_cra_H

#include "givaro/givconfig.h"
#include "givaro/giverror.h"

namespace Givaro {

	//! CRA
template<class Ring, class Domain, bool REDUCE = true>
struct ChineseRemainder {
    typedef typename Ring::Element   RingElement;
    typedef typename Domain::Element DomainElement;


        // Computes u = M^-1 in the domain
	// Then C_12 = (M^-1 mod D)M
    ChineseRemainder(const Ring& R, const RingElement& M, const Domain& D)
            : _domain(D)  {
        DomainElement u;
        _domain.invin( _domain.init(u, M) );
        _domain.convert(C_12, u);
        C_12 *= M;
    }



        // Computes res = A + ((e-A) mod D)*(M^-1 mod D)M
        // Then res mod M == A
        // And  res mod D == e
    RingElement & operator()( RingElement& res, const RingElement& A, const DomainElement& e) const {
        DomainElement smallA, smallM;
        _domain.init(smallA, A);
        _domain.init(smallM);
        _domain.sub(smallM, e, smallA);
        _domain.convert(res, smallM);
        res *= C_12;
        return res += A;
    }

private:
    Domain _domain;
    RingElement C_12;
};

//! CRA2.
//! JGD 05.12.2007: not required anymore ...

template<class Ring, class Domain>
struct ChineseRemainder<Ring, Domain, false>  {
    typedef typename Ring::Element   RingElement;
    typedef typename Domain::Element DomainElement;


        // Computes u = M^-1 in the domain
	// Then C_12 = (M^-1 mod D)M
    ChineseRemainder(const Ring& R, const RingElement& M, const Domain& D)
            : _domain(D) {
        DomainElement u;
        _domain.invin( _domain.init(u, M) );
        _domain.convert(C_12, u);
        C_12 *= M;
    }


        // Computes res = A + (e-A)*(M^-1 mod D)M
        // Then res mod M == A
        // And  res mod D == e
    RingElement & operator()( RingElement& res, const RingElement& A, const DomainElement& e) const {
        _domain.convert(res, e);
        res -= A;
        res *= C_12;
        return res += A;
    }

private:
    Domain _domain;
    RingElement C_12;
};

} // namespace Givaro
#endif // __GIVARO_cra_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
