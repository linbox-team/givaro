// ==========================================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// file: givadicqfq.h
// Time-stamp: <30 Nov 11 11:12:05 Jean-Guillaume.Dumas@imag.fr>
// date: 2007
// version:
// author: Jean-Guillaume.Dumas

/*! @file zpz/givtablelimits.h
 * @ingroup zpz
 * @brief  Zech extension.
 * Zech extension fitting a small enough memory space
 *   t-adic max sizes for BLAS based linear algebra over extension fields
 * @bib
 * - Dumas, Gautier, Pernet  <i>Finite field linear algebra subroutines.</i>
 *  ISSAC'02: Proceedings of the 2002 International Symposium on Symbolic
 * and Algebraic Computation, Lille, France pp 63--74.
 */

#ifndef __GIVARO_tablesize_MAX_H
#define __GIVARO_tablesize_MAX_H


#ifndef _GIVARO_FF_TABLE_MAX
// 2^23 ---> 2^23*4*3 = 100K
// #define FF_TABLE_MAX 8388608U
// 2^20 ---> 2s on 735MHz
//#define FF_TABLE_MAX 1048576U
// Now 2^21+1 seems OK
#define _GIVARO_FF_TABLE_MAX 2097153U
#endif

#ifndef _GIVARO_FF_MAXEXPONENT_
#define _GIVARO_FF_MAXEXPONENT_ 21
#endif


#include <iostream>
#include <vector>
#include "givaro/givprimes16.h"

#include <cmath>
#include <stddef.h>

namespace Givaro {

  // ---------------------------------------------  class
class AdicSize {
public:

    static size_t nmax53(const unsigned long P, const unsigned long e) {
        size_t i = 0;
        while (Primes16::ith(i) < P) ++i;
        return n_max_53[i][e-2];
    }

    static size_t qmax53(const unsigned long P, const unsigned long e) {
        size_t i = 0;
        while (Primes16::ith(i) < P) ++i;
        return qadic_53[i][e-2];
    }

    static size_t nmax64(const unsigned long P, const unsigned long e) {
        size_t i = 0;
        while (Primes16::ith(i) < P) ++i;
        return n_max_64[i][e-2];
    }

    size_t qmax64(const unsigned long P, const unsigned long e) {
        size_t i = 0;
        while (Primes16::ith(i) < P) ++i;
        return qadic_64[i][e-2];
    }

    static size_t twopmax53(const unsigned long P, const unsigned long e, const unsigned long nm) {
        double tmp = double(P-1);
        tmp *= double(P-1);
        tmp *= double(e);
        tmp *= double(nm);
        size_t k = (size_t)ilogb(tmp);
        return ( (53/(2*e-1))>k ? ++k : 0);
    }

    static size_t twopmax53(const unsigned long P, const unsigned long e) {
        size_t k = 53/(2*e-1);
        return ( std::pow((double)2.0,double(k))>double(e*(P-1)*(P-1)) ? k: 0);
    }

private:
    static const size_t n_max_53[][_GIVARO_FF_MAXEXPONENT_];
    static const size_t n_max_64[][_GIVARO_FF_MAXEXPONENT_];
    static const size_t qadic_53[][_GIVARO_FF_MAXEXPONENT_];
    static const size_t qadic_64[][_GIVARO_FF_MAXEXPONENT_];
};

} // namespace Givaro

#endif // __GIVARO_tablesize_MAX_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
