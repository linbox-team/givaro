#ifndef _CRA_H
#define _CRA_H
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givcra.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givcra.h,v 1.5 2005-02-02 19:08:29 pernet Exp $
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
            : _D(D)  {
      RingElement p;
      D.characteristic(p);
      P_12=  M * p ;
        DomainElement u;
        _D.invin( _D.init(u, M) );
        _D.convert(C_12, u);
        C_12 *= M;
    }



    RingElement & operator()( RingElement& res, const RingElement& A, const DomainElement& e) {
        _D.convert(res, e);
        res -= A;
        res *= C_12;
        res %= P_12;            // This is the reduction
        return res += A;
    }

private:
    Domain _D;
    RingElement C_12, P_12;
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
        _D.convert(C_12, u);
        C_12 *= M;
    }



    RingElement & operator()( RingElement& res, const RingElement& A, const DomainElement& e) {
        _D.convert(res, e);
        res -= A;
        res *= C_12;
        return res += A;
    }

private:
    Domain _D;
    RingElement C_12;
};


#endif
