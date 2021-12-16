// ==========================================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// file: givpowers.h
// Time-stamp: <15 Jun 15 16:36:32 Jean-Guillaume.Dumas@imag.fr>
// ==========================================================================

/** @file givpower.h
 * @ingroup system
 * @brief NO DOC
 */

#ifndef __GIVARO_power_H
#define __GIVARO_power_H

namespace Givaro {

    // -------------------------------------------------------------
    //! Integer log
    // -------------------------------------------------------------
    template<typename T>
    inline unsigned GIVINTLOG(const T& a)
    {
        unsigned l(0);
        for(T v(a); v >>= 1; ++l) {}
        return l;
    }



    // -------------------------------------------------------------
    //! Powering
    // -------------------------------------------------------------
    template<class TT, class UU>
    TT power(const TT n, const UU l)
    {
        if (l == 0) return 1 ;

        unsigned long p = (unsigned long) l ;
        short is_assg = 0 ;

        TT res = TT(1) ;
        TT puiss  = n ;

        while (p != 0) {
            if (p & 0x1) {
                if (is_assg)
                    res *= puiss ;
                else {
                    is_assg = 1 ;
                    res = puiss ;
                }
            }
            //       if ((p >>= 1) != 0) puiss = puiss * puiss ;
            if ((p >>= 1) != 0) puiss *= puiss ;

        }
        return res ;
    }

#if 0
#include <givaro/givinteger.h>
    template<> Integer power(const Integer n, const long l) { return pow(n,l); }
    template<> Integer power(const Integer n, const unsigned long l) { return pow(n,l); }
    template<> Integer power(const Integer n, const int l) { return pow(n,l); }
    template<> Integer power(const Integer n, const unsigned int l) { return pow(n,l); }
#endif

    //! dom_power
    template<class D, class TT>
    TT& dom_power(TT& res, const TT& n, uint64_t l, const D& F)
    {
        if (l == 0) return F.assign(res,F.one) ;

        uint64_t p = l;
        bool is_assg = false ;

        TT puiss; F.init(puiss); F.assign(puiss,n) ;
        F.assign(res,F.one) ;

        while (p != 0) {
            if (p & 0x1) {
                if (is_assg)
                    F.mulin(res,puiss) ;
                else {
                    is_assg = true ;
                    F.assign(res,puiss) ;
                }
            }
            if ((p >>= 1) != 0) { F.mulin(puiss,puiss) ; }
        }
        return res ;
    }


#if 0
#include <cmath>

    template<> double power<double>(const double a, const double e) {
        return pow(a,e);
    }
#endif

} // namespace Givaro

#endif // __GIVARO_power_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
