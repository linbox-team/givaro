#ifndef _CRA_H
#define _CRA_H
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givcra.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givcra.h,v 1.1 2004-12-17 14:42:59 jgdumas Exp $
// ==========================================================================
// Description:
//  Chinese Remainder Algorithm for 2 elements. 
//  For any number of moduli see givrns.h

#include "givaro/givconfig.h"
#include "givaro/giverror.h"


template<class Ring, class Domain, bool REDUCE>
struct ChineseRemainder;


template<>
template<class Ring, class Domain>
struct ChineseRemainder<Ring, Domain, false>  {
    typedef typename Ring::Element   RingElement;
    typedef typename Domain::Element DomainElement;


    ChineseRemainder(const Ring& R, const RingElement& M, const Domain& D) 
            : _D(D) {

        DomainElement u;
        _D.invin( _D.init(u, M) );
        _D.convert(_C12, u);
        _C12 *= M;
    }



    RingElement & operator()( RingElement& res, const RingElement& A, const DomainElement& e) {
        _D.convert(res, e);
        res -= A;
        res *= _C12;
        return res += A;
    }

private:
    Domain _D;
    RingElement _C12;

};

template<>
template<class Ring, class Domain>
struct ChineseRemainder<Ring, Domain, true>  {
    typedef typename Ring::Element   RingElement;
    typedef typename Domain::Element DomainElement;


    ChineseRemainder(const Ring& R, const RingElement& M, const Domain& D) 
            : _D(D), _P( M * (RingElement)D.characteristic() ) {

        DomainElement u;
        _D.invin( _D.init(u, M) );
        _D.convert(_C12, u);
        _C12 *= M;
    }



    RingElement & operator()( RingElement& res, const RingElement& A, const DomainElement& e) {
        _D.convert(res, e);
        res -= A;
        res *= _C12;
        res %= _P;
        return res += A;
    }

private:
    Domain _D;
    RingElement _C12, _P;

};

#endif
