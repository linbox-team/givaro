#ifndef _CRA_H
#define _CRA_H
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givcra.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givcra.h,v 1.10 2007-12-05 15:45:21 jgdumas Exp $
// ==========================================================================
// Description:
//  Chinese Remainder Algorithm for 2 Elements. 
//  For any number of moduli see givrns.h

#include "givaro/givconfig.h"
#include "givaro/giverror.h"


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

// JGD 05.12.2007: not required anymore ...
//template<>
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


#endif
