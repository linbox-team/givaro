#ifndef _CRA_H
#define _CRA_H
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givcra.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givcra.h,v 1.3 2005-01-06 17:10:50 jgdumas Exp $
// ==========================================================================
// Description:
//  Chinese Remainder Algorithm for 2 elements. 
//  For any number of moduli see givrns.h

#include "givaro/givconfig.h"
#include "givaro/giverror.h"


template<class Ring, class Domain, bool REDUCE = true>
struct ChineseRemainder {
    typedef typename Ring::Element   RingElement;
    typedef typename Domain::Element DomainElement;


    ChineseRemainder(const Ring& R, const RingElement& M, const Domain& D) 
            : _D(D), _P( M * (RingElement)D.characteristic() ) {

        DomainElement u;
        _D.invin( _D.init(u, M) );
        _D.convert(_C, u);
        _C *= M;
    }



    RingElement & operator()( RingElement& res, const RingElement& A, const DomainElement& e) {
        _D.convert(res, e);
        res -= A;
        res *= _C;
        res %= _P;            // This is the reduction
        return res += A;
    }

private:
    Domain _D;
    RingElement _C, _P;

};

template<>
template<class Ring, class Domain>
struct ChineseRemainder<Ring, Domain, false>  {
    typedef typename Ring::Element   RingElement;
    typedef typename Domain::Element DomainElement;


    ChineseRemainder(const Ring& R, const RingElement& M, const Domain& D) 
            : _D(D) {

        DomainElement u;
        _D.invin( _D.init(u, M) );
        _D.convert(_C, u);
        _C *= M;
    }



    RingElement & operator()( RingElement& res, const RingElement& A, const DomainElement& e) {
        _D.convert(res, e);
        res -= A;
        res *= _C;
        return res += A;
    }

private:
    Domain _D;
    RingElement _C;

};


#endif
