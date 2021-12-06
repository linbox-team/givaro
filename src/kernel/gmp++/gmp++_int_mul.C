// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int_mul.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: M. Samama, T. Gautier
// $Id: gmp++_int_mul.C,v 1.11 2011-01-20 08:19:15 jgdumas Exp $
// ==========================================================================
/** @file gmp++/gmp++_int_mul.C
 * muling stuff.
 */

#ifndef __GIVARO_gmpxx_gmpxx_int_mul_C
#define __GIVARO_gmpxx_gmpxx_int_mul_C

#ifndef __GIVARO_INLINE_ALL
#include "gmp++/gmp++.h"
#endif

namespace Givaro {

    //-------------------------------------------------- operator *
    Integer& Integer::mulin(Integer& res, const Integer& n)
    {
        if (isZero(n)) return res = Integer::zero;
        if (isZero(res)) return res;
        mpz_mul( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&n.gmp_rep );
        return res;
    }
    Integer& Integer::mulin(Integer& res, const int64_t n)
    {
        if (isZero(n)) return res = Integer::zero;
        if (isZero(res)) return res;
#if __GIVARO_SIZEOF_LONG < 8
        return mulin(res,Integer(n));
#else
        mpz_mul_si( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, n);
        return res;
#endif
    }
    Integer& Integer::mulin(Integer& res, const uint64_t n)
    {
        if (isZero(n)) return res = Integer::zero;
        if (isZero(res)) return res;
#if __GIVARO_SIZEOF_LONG < 8
        return mulin(res,Integer(n));
#else
        mpz_mul_ui( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&res.gmp_rep, n);
        return res;
#endif
    }

    Integer& Integer::mul(Integer& res, const Integer& n1, const Integer& n2)
    {
        if (isZero(n1)) return res = Integer::zero;
        if (isZero(n2)) return res = Integer::zero;
        mpz_mul( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&n1.gmp_rep, (mpz_srcptr)&n2.gmp_rep);
        return res;
    }
    Integer& Integer::mul(Integer& res, const Integer& n1, const int64_t n2)
    {
        if (isZero(n1)) return res = Integer::zero;
        if (isZero(n2)) return res = Integer::zero;
#if __GIVARO_SIZEOF_LONG < 8
        return mul(res,n1,Integer(n2));
#else
        mpz_mul_si( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&n1.gmp_rep, n2);
        return res;
#endif
    }
    Integer& Integer::mul(Integer& res, const Integer& n1, const uint64_t n2)
    {
        if (isZero(n1)) return res = Integer::zero;
        if (isZero(n2)) return res = Integer::zero;
#if __GIVARO_SIZEOF_LONG < 8
        return mul(res,n1,Integer(n2));
#else
        mpz_mul_ui( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&n1.gmp_rep, n2);
        return res;
#endif
    }

    Integer& Integer::axpy(Integer& res, const Integer& a, const Integer& x, const Integer& b)
    {
        if (&res == &b) return Integer::axpyin(res,a,x);
        if (isZero(a) || isZero(x)) return res = b;
        mpz_mul( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&a.gmp_rep, (mpz_srcptr)&x.gmp_rep);
        mpz_add( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&b.gmp_rep);
        return res;
    }

    Integer& Integer::axpy(Integer& res, const Integer& a, const uint64_t x, const Integer& b)
    {
        if (&res == &b) return Integer::axpyin(res,a,x);
        if (isZero(a) || isZero(x)) return res = b;
#if __GIVARO_SIZEOF_LONG < 8
        return axpy(res,a,Integer(x),b);
#else
        mpz_mul_ui( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&a.gmp_rep, x);
        mpz_add( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&b.gmp_rep);
        return res;
#endif
    }


    Integer& Integer::axpyin(Integer& res, const Integer& a, const Integer& x)
    {
        if (isZero(a) || isZero(x)) return res;
        mpz_addmul( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&a.gmp_rep, (mpz_srcptr)&x.gmp_rep);
        return res;
    }

    Integer& Integer::axpyin(Integer& res, const Integer& a, const uint64_t x)
    {
        if (isZero(a) || isZero(x)) return res;
#if __GIVARO_SIZEOF_LONG < 8
        return axpyin(res,a,Integer(x));
#else
        mpz_addmul_ui( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&a.gmp_rep, x);
        return res;
#endif
    }


    Integer& Integer::maxpy(Integer& res, const Integer& a, const Integer& x, const Integer& b)
    {
        if (isZero(a) || isZero(x)) return res=b;
        if (&res == &b) return Integer::maxpyin(res,a,x);
        mpz_mul( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&a.gmp_rep, (mpz_srcptr)&x.gmp_rep);
        mpz_sub( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&b.gmp_rep, (mpz_ptr)&res.gmp_rep);
        return res;
    }
    Integer& Integer::maxpy(Integer& res, const Integer& a, const uint64_t x, const Integer& b)
    {
        if (isZero(a) || isZero(x)) return res=b;
        if (&res == &b) return Integer::maxpyin(res,a,x);
#if __GIVARO_SIZEOF_LONG < 8
        return maxpy(res,a,Integer(x),b);
#else
        mpz_mul_ui( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&a.gmp_rep, x);
        mpz_sub( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&b.gmp_rep, (mpz_ptr)&res.gmp_rep);
        return res;
#endif
    }

    Integer& Integer::axmy(Integer& res, const Integer& a, const Integer& x, const Integer& b)
    {
        if (&res == &b) return Integer::axmyin(res,a,x);
        if (isZero(a) || isZero(x)) return Integer::neg(res,b);
        mpz_mul( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&a.gmp_rep, (mpz_srcptr)&x.gmp_rep);
        mpz_sub( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&b.gmp_rep);
        return res;
    }

    Integer& Integer::axmy(Integer& res, const Integer& a, const uint64_t x, const Integer& b)
    {
        if (&res == &b) return Integer::axmyin(res,a,x);
        if (isZero(a) || isZero(x)) return Integer::neg(res,b);
#if __GIVARO_SIZEOF_LONG < 8
        return axmy(res,a,Integer(x),b);
#else
        mpz_mul_ui( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&a.gmp_rep, x);
        mpz_sub( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&b.gmp_rep);
        return res;
#endif
    }

    Integer& Integer::axmyin(Integer& res, const Integer& a, const Integer& x)
    {
        Integer::maxpyin(res,a,x);
        Integer::negin(res);
        return res ;
    }

    Integer& Integer::axmyin(Integer& res, const Integer& a, const uint64_t x)
    {
        Integer::maxpyin(res,a,x);
        Integer::negin(res);
        return res ;
    }

    Integer& Integer::maxpyin(Integer& res, const Integer& a, const Integer& x)
    {
        if (isZero(a) || isZero(x)) return res;
        mpz_submul( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&a.gmp_rep, (mpz_srcptr)&x.gmp_rep);
        return res;
    }

    Integer& Integer::maxpyin(Integer& res, const Integer& a, const uint64_t x)
    {
        if (isZero(a) || isZero(x)) return res;
#if __GIVARO_SIZEOF_LONG < 8
        return maxpyin(res,a,Integer(x));
#else
        mpz_submul_ui( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&a.gmp_rep, x);
        return res;
#endif
    }

    Integer& Integer::operator *= (const Integer& n)
    {
        if (isZero(n)) return *this = Integer::zero;
        if (isZero(*this)) return *this;
        //   Rep (res.gmp_rep)( MAX(SZ_REP(n.gmp_rep),SZ_REP(gmp_rep)) );
        Integer res;
        mpz_mul( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, (mpz_srcptr)&n.gmp_rep) ;
        return *this = res;
    }

    Integer& Integer::operator *= (const uint64_t l)
    {
        if (l==0) return *this = Integer::zero;
        if (isZero(*this)) return *this;
#if __GIVARO_SIZEOF_LONG < 8
        return mulin(*this,Integer(l));
#else
        mpz_mul_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, l);
        return *this;
#endif
    }

    Integer& Integer::operator *= (const int64_t l)
    {
        if (l==0) return *this =Integer::zero;
        if (isZero(*this)) return *this;
#if __GIVARO_SIZEOF_LONG < 8
        return mulin(*this,Integer(l));
#else
        mpz_mul_si( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, l);
        return *this;
#endif
    }


    Integer Integer::operator * (const Integer& n) const
    {
        if (isZero(n)) return Integer::zero;
        if (isZero(*this)) return Integer::zero;
        Integer res;
        mpz_mul( (mpz_ptr)&(res.gmp_rep), (mpz_srcptr)&gmp_rep, (mpz_srcptr)&n.gmp_rep) ;
        return res;
    }

    Integer Integer::operator * (const uint64_t l) const
    {
        if (l==0) return Integer::zero;
        if (isZero(*this)) return Integer::zero;
#if __GIVARO_SIZEOF_LONG < 8
        return (*this) * Integer(l);
#else
        Integer res;
        mpz_mul_ui( (mpz_ptr)&(res.gmp_rep), (mpz_srcptr)&gmp_rep, l);
        return res;
#endif
    }

    Integer Integer::operator * (const int64_t l) const
    {
        if (l==0) return Integer::zero;
        if (isZero(*this)) return Integer::zero;
#if __GIVARO_SIZEOF_LONG < 8
        return (*this) * Integer(l);
#else
        Integer res;
        mpz_mul_si( (mpz_ptr)&(res.gmp_rep), (mpz_srcptr)&gmp_rep, l);
        return res;
#endif
    }

    // -- operator *
    Integer operator * (const int32_t l, const Integer& n)
    {
        return n * (int64_t)l;
    }
    Integer operator * (const uint32_t l, const Integer& n)
    {
        return n * (uint64_t)l;
    }
    Integer operator * (const int64_t l, const Integer& n)
    {
        return n * l;
    }
    Integer operator * (const uint64_t l, const Integer& n)
    {
        return n * l;
    }



}
#endif // __GIVARO_gmpxx_gmpxx_int_mul_C
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
