// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int_div.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: M. Samama, T. Gautier
// $Id: gmp++_int_div.C,v 1.9 2011-01-06 18:02:37 briceboyer Exp $
// ==========================================================================
/** @file gmp++/gmp++_int_div.C
 * diving stuff.
 */

#ifndef __GIVARO_gmpxx_gmpxx_int_div_C
#define __GIVARO_gmpxx_gmpxx_int_div_C

#ifndef __GIVARO_INLINE_ALL
#include "gmp++/gmp++.h"
#endif
#include <cstdlib>

namespace Givaro {

    //-------------------------------------------------- operator /
    Integer& Integer::divin(Integer& res, const Integer& n)
    {
        //  if (isZero(n)) {
        //    GivMathDivZero("[Integer::/]: division by zero");
        //  }
        if (isZero(res)) return res;
        mpz_tdiv_q( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&n.gmp_rep );
        return res;
    }

    Integer& Integer::divin(Integer& res, const int64_t n)
    {

        if (isZero(res)) return res;
#if __GIVARO_SIZEOF_LONG < 8
        return divin(res,Integer(n));
#else
        int32_t sgn = Givaro::sign(n);
        mpz_tdiv_q_ui( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, std::abs(n));
        if (sgn <0) return res = -res;
        return res;
#endif
    }

    Integer& Integer::divin(Integer& res, const uint64_t n)
    {
        if (isZero(res)) return res;
#if __GIVARO_SIZEOF_LONG < 8
        return divin(res,Integer(n));
#else
        mpz_tdiv_q_ui( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&res.gmp_rep, n);
        return res;
#endif
    }

    Integer& Integer::div(Integer& res, const Integer& n1, const Integer& n2)
    {
        if (isZero(n1)) return res = Integer::zero;
        //  if (isZero(n2)) {
        //    GivMathDivZero("[Integer::/]: division by zero");
        //  }
        mpz_tdiv_q( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&n1.gmp_rep, (mpz_srcptr)&n2.gmp_rep);
        return res;
    }

    Integer& Integer::div(Integer& res, const Integer& n1, const int64_t n2)
    {
        if (isZero(n1)) return res = Integer::zero;
#if __GIVARO_SIZEOF_LONG < 8
        return div(res,n1,Integer(n2));
#else
        int32_t sgn = Givaro::sign(n2);
        mpz_tdiv_q_ui( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&n1.gmp_rep, std::abs(n2));
        if (sgn <0) return res = -res;
        return res;
#endif
    }

    Integer& Integer::div(Integer& res, const Integer& n1, const int32_t n2)
    {
        return div(res,n1,int64_t(n2));
    }

    Integer& Integer::div(Integer& res, const Integer& n1, const uint64_t n2)
    {
        if (isZero(n1)) return res = Integer::zero;
#if __GIVARO_SIZEOF_LONG < 8
        return div(res,n1,Integer(n2));
#else
        mpz_tdiv_q_ui( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&n1.gmp_rep, n2);
        return res;
#endif
    }

    Integer& Integer::divexact  (Integer& q, const Integer& n1, const Integer& n2)
    {
        if (isZero(n1)) return q = Integer::zero;
        //  if (isZero(n2)) {
        //    GivMathDivZero("[Integer::/]: division by zero");
        //  }
        mpz_divexact( (mpz_ptr)&(q.gmp_rep),
                      (mpz_srcptr)&(n1.gmp_rep), (mpz_srcptr)&(n2.gmp_rep)) ;
        return q;
    }

    Integer& Integer::divexact  (Integer& q, const Integer& n1, const uint64_t & n2)
    {
        if (isZero(n1)) return q = Integer::zero;
#if __GIVARO_SIZEOF_LONG < 8
        return divexact(q,n1,Integer(n2));
#else
        mpz_divexact_ui( (mpz_ptr)&(q.gmp_rep),
                         (mpz_srcptr)&(n1.gmp_rep), (n2)) ;
        return q;
#endif
    }

    Integer& Integer::divexact  (Integer& q, const Integer& n1, const int64_t& n2)
    {
        if (isZero(n1)) return q = Integer::zero;
#if __GIVARO_SIZEOF_LONG < 8
        return divexact(q,n1,Integer(n2));
#else
        mpz_divexact_ui( (mpz_ptr)&(q.gmp_rep),
                         (mpz_srcptr)&(n1.gmp_rep), std::abs(n2)) ;
        if (n2<0)
            negin(q);
        return q;
#endif
    }


    Integer  Integer::divexact  (const Integer& n1, const Integer& n2)
    {
        if (isZero(n1)) return Integer::zero;
        //  if (isZero(n2)) {
        //    GivMathDivZero("[Integer::/]: division by zero");
        //  }
        Integer q;
        mpz_divexact( (mpz_ptr)&(q.gmp_rep),
                      (mpz_srcptr)&(n1.gmp_rep), (mpz_srcptr)&(n2.gmp_rep)) ;
        return q;
    }

    Integer  Integer::divexact  (const Integer& n1, const uint64_t& n2)
    {
        if (isZero(n1)) return Integer::zero;
#if __GIVARO_SIZEOF_LONG < 8
        return divexact(n1,Integer(n2));
#else
        Integer q;
        mpz_divexact_ui( (mpz_ptr)&(q.gmp_rep),
                         (mpz_srcptr)&(n1.gmp_rep), (n2)) ;
        return q;
#endif
    }

    Integer  Integer::divexact  (const Integer& n1, const int64_t& n2)
    {
        if (isZero(n1)) return Integer::zero;
#if __GIVARO_SIZEOF_LONG < 8
        return divexact(n1,Integer(n2));
#else
        Integer q;
        mpz_divexact_ui( (mpz_ptr)&(q.gmp_rep),
                         (mpz_srcptr)&(n1.gmp_rep), std::abs(n2)) ;
        if (n2<0) negin(q);
        return q;
#endif
    }

    Integer& Integer::operator /= (const Integer& n)
    {
        //  if (isZero(n)) {
        //    GivMathDivZero("[Integer::/]: division by zero");
        //  }
        if (isZero(*this)) return *this;
        mpz_tdiv_q( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, (mpz_srcptr)&n.gmp_rep) ;
        return *this;
    }

    Integer& Integer::operator /= (const uint64_t l)
    {
        //  if (l ==0) {
        //    GivMathDivZero("[Integer::/]: division by zero");
        //  }
        if (isZero(*this)) return *this;
#if __GIVARO_SIZEOF_LONG < 8
        return *this /= Integer(l);
#else
        mpz_tdiv_q_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, l);
        return *this;
#endif
    }

    Integer& Integer::operator /= (const int64_t l)
    {
        //  if (l ==0) {
        //    GivMathDivZero("[Integer::/]: division by zero");
        //  }
        if (isZero(*this)) return *this;
#if __GIVARO_SIZEOF_LONG < 8
        return *this /= Integer(l);
#else
        int32_t sgn = Givaro::sign(l);
        mpz_tdiv_q_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, std::abs(l));
        if (sgn <0) mpz_neg( (mpz_ptr)&gmp_rep, (mpz_ptr)&(gmp_rep));
        return *this;
#endif
    }


    Integer Integer::operator / (const Integer& n) const
    {
        //  if (isZero(n)) {
        //    GivMathDivZero("[Integer::/]: division by zero");
        //  }
        if (isZero(*this)) return Integer::zero;
        Integer res;
        mpz_tdiv_q( (mpz_ptr)&(res.gmp_rep), (mpz_srcptr)&gmp_rep, (mpz_srcptr)&n.gmp_rep) ;
        return res;
    }

    Integer Integer::operator / (const uint64_t l) const
    {
        //  if (l ==0) {
        //    GivMathDivZero("[Integer::/]: division by zero");
        //  }
        if (isZero(*this)) return Integer::zero;
#if __GIVARO_SIZEOF_LONG < 8
        return *this / Integer(l);
#else
        Integer res;
        mpz_tdiv_q_ui( (mpz_ptr)&(res.gmp_rep), (mpz_srcptr)&gmp_rep, l);

        return res;
#endif
    }

    Integer Integer::operator / (const int64_t l) const
    {
        //  if (l ==0) {
        //    GivMathDivZero("[Integer::/]: division by zero");
        //  }
        if (isZero(*this)) return Integer::zero;
#if __GIVARO_SIZEOF_LONG < 8
        return *this / Integer(l);
#else
        Integer res;
        int32_t sgn = Givaro::sign(l);
        mpz_tdiv_q_ui( (mpz_ptr)&(res.gmp_rep), (mpz_srcptr)&gmp_rep, std::abs(l));
        if (sgn <0) return negin(res);
        return res;
#endif
    }

    // Euclidian division
    Integer& Integer::divmod(Integer& q, Integer& r, const Integer &a, const Integer &b)
    {
        mpz_tdiv_qr( (mpz_ptr)&(q.gmp_rep), (mpz_ptr)&(r.gmp_rep),
                     (mpz_srcptr)&(a.gmp_rep), (mpz_srcptr)&(b.gmp_rep));

        /* If r is negative (happen if a is negative, as sign(r) = sign(a)),
         * we need to modify q and r to have a positive r.
         */
        if (r < 0)
        {
            if (b > 0)
            {
                subin (q, (uint64_t)1) ;
                r += b;
            }
            else /* b is negative */
            {
                addin (q, (uint64_t)1) ;
                r -= b;
            }
        }

        return q;
    }

    Integer& Integer::divmod(Integer& q, int64_t & r, const Integer& a, const int64_t b)
    {
#if __GIVARO_SIZEOF_LONG < 8
        Integer res;
        divmod(q,res,a,Integer(b));
        r = (int64_t)res;
        return q;
#else
        r = (int64_t)mpz_tdiv_q_ui( (mpz_ptr)&(q.gmp_rep),
                                    (mpz_srcptr)&(a.gmp_rep), std::abs(b));

        if (a<0 && r) {
            // :GMPUintTDiv
            subin(q,(int64_t)1) ;
            r = b - r ;
        }

        return q;
#endif
    }

    Integer& Integer::divmod(Integer& q, uint64_t & r, const Integer& a, const uint64_t b)
    {
#if __GIVARO_SIZEOF_LONG < 8
        Integer res;
        divmod(q,res,a,Integer(b));
        r = (uint64_t)res; // divmod already corrects when a<0
        return q;
#else
        r = mpz_tdiv_q_ui( (mpz_ptr)&(q.gmp_rep),
                           (mpz_srcptr)&(a.gmp_rep), b);

        if (a<0 && r) {
            subin(q,(int64_t)1) ;
            // :GMPUintTDiv The GMP documentation specifies that:
            // 'For the ui variants (...) tdiv and cdiv the remainder can be negative,
            // so for those the return value is the absolute value of the remainder.'
            // Thus we have to correct this with (-q1-1),d-r1
            r = b - r;
        }

        return q;
#endif
    }

    Integer& Integer::ceil(Integer& q, const Integer & n, const Integer & d)
    {
        mpz_cdiv_q( (mpz_ptr)&(q.gmp_rep),
                    (mpz_srcptr)&(n.gmp_rep),
                    (mpz_srcptr)&(d.gmp_rep));
        return q ;
    }

    Integer& Integer::floor(Integer& q, const Integer & n, const Integer & d)
    {
        mpz_fdiv_q( (mpz_ptr)&(q.gmp_rep),
                    (mpz_srcptr)&(n.gmp_rep),
                    (mpz_srcptr)&(d.gmp_rep));
        return q ;
    }

    Integer& Integer::trunc(Integer& q, const Integer & n, const Integer & d)
    {
        mpz_tdiv_q( (mpz_ptr)&(q.gmp_rep),
                    (mpz_srcptr)&(n.gmp_rep),
                    (mpz_srcptr)&(d.gmp_rep));
        return q ;
    }

    Integer Integer::ceil( const Integer & n, const Integer & d)
    {
        Integer q;
        mpz_cdiv_q( (mpz_ptr)&(q.gmp_rep),
                    (mpz_srcptr)&(n.gmp_rep),
                    (mpz_srcptr)&(d.gmp_rep));
        return q ;
    }

    Integer Integer::floor(const Integer & n, const Integer & d)
    {
        Integer q;
        mpz_fdiv_q( (mpz_ptr)&(q.gmp_rep),
                    (mpz_srcptr)&(n.gmp_rep),
                    (mpz_srcptr)&(d.gmp_rep));
        return q ;
    }

    Integer Integer::trunc(const Integer & n, const Integer & d)
    {
        Integer q;
        mpz_tdiv_q( (mpz_ptr)&(q.gmp_rep),
                    (mpz_srcptr)&(n.gmp_rep),
                    (mpz_srcptr)&(d.gmp_rep));
        return q ;
    }

    Integer& Integer::trem(Integer& r, const Integer &n , const Integer & d)
    {
        mpz_tdiv_r((mpz_ptr)&(r.gmp_rep),
                   (mpz_srcptr)&(n.gmp_rep),
                   (mpz_srcptr)&(d.gmp_rep));
        return r;
    }

    Integer& Integer::crem(Integer& r, const Integer &n , const Integer & d)
    {
        mpz_cdiv_r((mpz_ptr)&(r.gmp_rep),
                   (mpz_srcptr)&(n.gmp_rep),
                   (mpz_srcptr)&(d.gmp_rep));
        return r;
    }

    Integer& Integer::frem(Integer& r, const Integer &n , const Integer & d)
    {
        mpz_fdiv_r((mpz_ptr)&(r.gmp_rep),
                   (mpz_srcptr)&(n.gmp_rep),
                   (mpz_srcptr)&(d.gmp_rep));
        return r;
    }

    Integer& Integer::trem(Integer& r, const Integer &n , const uint64_t& d)
    {
#if __GIVARO_SIZEOF_LONG < 8
        return trem(r,n,Integer(d));
#else
        mpz_tdiv_r_ui((mpz_ptr)&(r.gmp_rep),
                      (mpz_srcptr)&(n.gmp_rep),
                      (d));
        return r;
#endif
    }

    Integer& Integer::crem(Integer& r, const Integer &n , const uint64_t & d)
    {
#if __GIVARO_SIZEOF_LONG < 8
        return crem(r,n,Integer(d));
#else
        mpz_cdiv_r_ui((mpz_ptr)&(r.gmp_rep),
                      (mpz_srcptr)&(n.gmp_rep),
                      (d));
        return r;
#endif
    }

    Integer& Integer::frem(Integer& r, const Integer &n , const uint64_t & d)
    {
#if __GIVARO_SIZEOF_LONG < 8
        return frem(r,n,Integer(d));
#else
        mpz_fdiv_r_ui((mpz_ptr)&(r.gmp_rep),
                      (mpz_srcptr)&(n.gmp_rep),
                      d);
        return r;
#endif
    }

    uint64_t Integer::trem(const Integer &n , const uint64_t& d)
    {
#if __GIVARO_SIZEOF_LONG < 8
        return (uint64_t)trem(n,Integer(d));
#else
        return mpz_cdiv_ui( (mpz_srcptr)&(n.gmp_rep),
                            (d));
#endif
    }

    uint64_t Integer::crem(const Integer &n , const uint64_t & d)
    {
#if __GIVARO_SIZEOF_LONG < 8
        return (uint64_t)crem(n,Integer(d));
#else
        return mpz_tdiv_ui( (mpz_srcptr)&(n.gmp_rep),
                            (d));
#endif
    }

    uint64_t Integer::frem(const Integer &n , const uint64_t & d)
    {
#if __GIVARO_SIZEOF_LONG < 8
        return (uint64_t)frem(n,Integer(d));
#else
        return mpz_fdiv_ui( (mpz_srcptr)&(n.gmp_rep),
                            d);
#endif
    }

    // -- operator /
    Integer operator / (const int32_t l, const Integer& n)
    {
        return Integer(l)/n;
    }
    Integer operator / (const int64_t l, const Integer& n)
    {
        return Integer(l)/n;
    }
    // -- operator /
    Integer operator / (const uint32_t l, const Integer& n)
    {
        return Integer(l)/n;
    }
    Integer operator / (const uint64_t l, const Integer& n)
    {
        return Integer(l)/n;
    }



}

#endif // __GIVARO_gmpxx_gmpxx_int_div_C

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
