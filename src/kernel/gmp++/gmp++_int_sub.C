// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int_sub.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: M. Samama, T. Gautier
// $Id: gmp++_int_sub.C,v 1.4 2009-09-17 14:28:22 jgdumas Exp $
// ==========================================================================

/** @file gmp++/gmp++_int_sub.C
 * subing stuff.
 */

#ifndef __GIVARO_gmpxx_gmpxx_int_sub_C
#define __GIVARO_gmpxx_gmpxx_int_sub_C

#ifndef __GIVARO_INLINE_ALL
#include "gmp++/gmp++.h"
#endif

namespace Givaro {

    //-------------------------------------------------- operator -
    Integer& Integer::subin(Integer& res, const Integer& n)
    {
        if (isZero(n)) return res;
        if (isZero(res)) return res = - n;
        mpz_sub( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&n.gmp_rep );
        return res;
    }
    Integer& Integer::subin(Integer& res, const int64_t n)
    {
        if (isZero(n)) return res;
        if (isZero(res)) return negin(res = n);
#if __GIVARO_SIZEOF_LONG < 8
        return subin(res,Integer(n));
#else
        int32_t sgn = Givaro::sign(n);
        if (sgn >0) mpz_sub_ui( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, n);
        else mpz_add_ui((mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, -static_cast<uint64_t>(n));
        return res;
#endif
    }
    Integer& Integer::subin(Integer& res, const uint64_t n)
    {
        if (isZero(n)) return res;
        if (isZero(res)) return negin(res = n);
#if __GIVARO_SIZEOF_LONG < 8
        return subin(res,Integer(n));
#else
        mpz_sub_ui( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&res.gmp_rep, n);
        return res;
#endif
    }

    Integer& Integer::sub(Integer& res, const Integer& n1, const Integer& n2)
    {
        if (isZero(n1)) return res = - n2;
        if (isZero(n2)) return res = n1;
        mpz_sub( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&n1.gmp_rep, (mpz_srcptr)&n2.gmp_rep);
        return res;
    }
    Integer& Integer::sub(Integer& res, const Integer& n1, const int64_t n2)
    {
        if (isZero(n1)) return negin(res = n2);
        if (isZero(n2)) return res = n1;
#if __GIVARO_SIZEOF_LONG < 8
        return sub(res,n1,Integer(n2));
#else
        int32_t sgn = Givaro::sign(n2);
        if (sgn >0) mpz_sub_ui( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&n1.gmp_rep, n2);
        else mpz_add_ui((mpz_ptr)&res.gmp_rep, (mpz_srcptr)&n1.gmp_rep, -static_cast<uint64_t>(n2));
        return res;
#endif
    }
    Integer& Integer::sub(Integer& res, const Integer& n1, const uint64_t n2)
    {
        if (isZero(n1)) return negin(res = n2);
        if (isZero(n2)) return res = n1;
#if __GIVARO_SIZEOF_LONG < 8
        return sub(res,n1,Integer(n2));
#else
        mpz_sub_ui( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&n1.gmp_rep, n2);
        return res;
#endif
    }

    Integer& Integer::neg(Integer& res, const Integer& n)
    {
        mpz_neg( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&n.gmp_rep);
        return res;
    }

    Integer& Integer::negin(Integer& res)
    {
        mpz_neg( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep);
        return res;
    }


    Integer& Integer::operator -= (const Integer& n)
    {
        if (isZero(n)) return *this;
        if (isZero(*this)) return logcpy(-n);
        //   Rep (res.gmp_rep)( MAX(SZ_REP(n.gmp_rep),SZ_REP(gmp_rep)) );
        mpz_sub( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, (mpz_srcptr)&n.gmp_rep) ;
        return *this;
    }

    Integer& Integer::operator -= (const uint64_t l)
    {
        if (l==0) return *this;
        if (isZero(*this)) return logcpy(-Integer(l));
#if __GIVARO_SIZEOF_LONG < 8
        return (*this) -= Integer(l);
#else
        //   Rep (res.gmp_rep)( MAX(SZ_REP(gmp_rep),1) );
        mpz_sub_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, l);
        return *this;
#endif
    }

    Integer& Integer::operator -= (const int64_t l)
    {
        if (l==0) return *this;
        if (isZero(*this)) return logcpy(-Integer(l));
        //   Rep (res.gmp_rep)( MAX(SZ_REP(gmp_rep),1) );
#if __GIVARO_SIZEOF_LONG < 8
        return (*this) -= Integer(l);
#else
        int32_t sgn = Givaro::sign(l);
        if (sgn >0) mpz_sub_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, l);
        else mpz_add_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, -static_cast<uint64_t>(l));
        return *this;
#endif
    }


    Integer Integer::operator - (const Integer& n) const
    {
        if (isZero(n)) return *this;
        if (isZero(*this)) return -n;
        //   Rep (res.gmp_rep)( MAX(SZ_REP(n.gmp_rep),SZ_REP(gmp_rep)) );
        Integer res;
        mpz_sub( (mpz_ptr)&(res.gmp_rep), (mpz_srcptr)&gmp_rep, (mpz_srcptr)&n.gmp_rep) ;
        return res;
    }

    Integer Integer::operator - (const uint64_t l) const
    {
        if (l==0) return *this;
        if (isZero(*this)) return -Integer(l);
        //   Rep (res.gmp_rep)( MAX(SZ_REP(gmp_rep),1) );
#if __GIVARO_SIZEOF_LONG < 8
        return (*this) - Integer(l);
#else
        Integer res;
        mpz_sub_ui( (mpz_ptr)&(res.gmp_rep), (mpz_srcptr)&gmp_rep, l);
        return res;
#endif
    }

    Integer Integer::operator - (const int64_t l) const
    {
        if (l==0) return *this;
        if (isZero(*this)) return -Integer(l);
        //   Rep (res.gmp_rep)( MAX(SZ_REP(gmp_rep),1) );
#if __GIVARO_SIZEOF_LONG < 8
        return (*this) - Integer(l);
#else
        Integer res;
        int32_t sgn = Givaro::sign(l);
        if (sgn >0) mpz_sub_ui( (mpz_ptr)&(res.gmp_rep), (mpz_srcptr)&gmp_rep, l);
        else mpz_add_ui( (mpz_ptr)&(res.gmp_rep), (mpz_srcptr)&gmp_rep, -static_cast<uint64_t>(l));
        return res;
#endif
    }

    // -- operator -
    Integer operator - (const int32_t l, const Integer& n)
    {
        return -(n - (int64_t)l);
    }
    Integer operator - (const uint32_t l, const Integer& n)
    {
        return -(n - (uint64_t)l);
    }
    Integer operator - (const int64_t l, const Integer& n)
    {
        return -(n - l);
    }
    Integer operator - (const uint64_t l, const Integer& n)
    {
        return -(n - l);
    }


}

#endif // __GIVARO_gmpxx_gmpxx_int_sub_C

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
