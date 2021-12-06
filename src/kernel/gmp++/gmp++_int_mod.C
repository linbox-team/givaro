// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int_mod.C,v $
// Copyright(c)'1994-2010 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: M. Samama, T. Gautier
// Modified: JG. Dumas, BB.
// $Id: gmp++_int_mod.C,v 1.17 2010-12-22 13:47:45 jgdumas Exp $
// ==========================================================================
/** @file gmp++/gmp++_int_mod.C
 * moding stuff.
 */

#ifndef __GIVARO_gmpxx_gmpxx_int_mod_C
#define __GIVARO_gmpxx_gmpxx_int_mod_C

#ifndef __GIVARO_INLINE_ALL
#include "gmp++/gmp++.h"
#endif
#include <cstdlib>

#ifndef GIVABS
#define GIVABS(a) ((a)>0?(a):-(a))
#endif

namespace Givaro {
    //-------------------------------------------------- operator /
    Integer& Integer::modin(Integer& res, const Integer& n)
    {
        if (isZero(res)) return res;
        //   mpz_tdiv_r( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n.gmp_rep );
        mpz_mod( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&n.gmp_rep );
        return res;
    }
    Integer& Integer::modin(Integer& res, const uint64_t n)
    {
        if (isZero(res)) return res;
#if __GIVARO_SIZEOF_LONG < 8
        return modin(res,Integer(n));
#else
        mpz_mod_ui( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&res.gmp_rep, n);
        return res;
#endif
    }
    Integer& Integer::modin(Integer& res, const int64_t n)
    {
        if (isZero(res)) return res;
#if __GIVARO_SIZEOF_LONG < 8
        return modin(res,Integer(n));
#else
        if (n>0)
            mpz_mod_ui( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&res.gmp_rep, n);
        else
            mpz_mod_ui( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&res.gmp_rep, -n);
        return res;
#endif
    }

    Integer& Integer::mod(Integer& res, const Integer& n1, const Integer& n2)
    {
        if (isZero(n1)) return res = Integer::zero;
        mpz_mod( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&n1.gmp_rep, (mpz_srcptr)&n2.gmp_rep);
        assert(!(res<0) && (res<abs(n2)));
        return res;
    }
    Integer& Integer::mod(Integer& res, const Integer& n1, const int64_t n2)
    {
        if (isZero(n1)) return res = Integer::zero;
#if __GIVARO_SIZEOF_LONG < 8
        return mod(res,n1,Integer(n2));
#else
        if (n2>0)
            mpz_mod_ui( (mpz_ptr)&(res.gmp_rep), (mpz_srcptr)&n1.gmp_rep, n2);
        else
            mpz_mod_ui( (mpz_ptr)&(res.gmp_rep), (mpz_srcptr)&n1.gmp_rep, -n2);

        assert(!(res<0) && (res<std::abs(n2)));
        return res;
#endif
    }
    Integer& Integer::mod(Integer& res, const Integer& n1, const uint64_t n2)
    {
        if (isZero(n1)) return res = Integer::zero;
#if __GIVARO_SIZEOF_LONG < 8
        return mod(res,n1,Integer(n2));
#else
        mpz_mod_ui( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&n1.gmp_rep, n2);
        assert(!(res<0) && (res<n2));
        return res;
#endif
    }


    Integer& Integer::operator %= (const Integer& n)
    {
        if (isZero(*this)) return *this;
        mpz_tdiv_r( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, (mpz_srcptr)&n.gmp_rep) ;
        return *this;
    }

    Integer& Integer::operator %= (const uint64_t l)
    {
        if (isZero(*this)) return *this;
#if __GIVARO_SIZEOF_LONG < 8
        return *this %= Integer(l);
#else
#  ifdef __GIVARO_DEBUG
        int32_t sgn_this = (*this>0)?1:-1;
#  endif
        mpz_tdiv_r_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, l);

#  ifdef __GIVARO_DEBUG
        assert((*this<GIVABS(l)) && (*this> -GIVABS(l)) && (sgn_this*(*this).priv_sign()>=0)) ;
#  endif
        return *this;
#endif
    }

    Integer& Integer::operator %= (const int64_t l)
    {
#if __GIVARO_SIZEOF_LONG < 8
        return *this %= Integer(l);
#else
        if (isZero(*this)) return *this;
#  ifdef __GIVARO_DEBUG
        int32_t sgn_this = (*this>0)?1:-1;
#  endif
        int32_t sgn = Givaro::sign(l);
        mpz_tdiv_r_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, std::abs(l));
        if (sgn <0) mpz_neg( (mpz_ptr)&gmp_rep, (mpz_ptr)&(gmp_rep) );

#  ifdef __GIVARO_DEBUG
        assert((*this<GIVABS(l)) && (*this> -GIVABS(l)) && (sgn_this*(*this).priv_sign()>=0)) ;
#  endif
        return *this;
#endif
    }

    Integer Integer::operator % (const Integer& n) const
    {
        if (isZero(*this)) return Integer::zero;
        Integer res;
        mpz_tdiv_r( (mpz_ptr)&(res.gmp_rep), (mpz_srcptr)&gmp_rep, (mpz_srcptr)&n.gmp_rep) ;

        // std::cout << res << ',' << n << ',' << *this << std::endl;
        assert((res<GIVABS(n)) && (res> -GIVABS(n)) && (res.priv_sign()*(*this).priv_sign()>=0)) ;
        return res;
    }

    int64_t Integer::operator % (const uint64_t l) const
    {
        if (isZero(*this)) return 0U;
#if __GIVARO_SIZEOF_LONG < 8
        return (int64_t) ((*this) % Integer(l));
#else
        bool isneg = (*this)<0 ;
        //CONDITION: mpz_tdiv_ui does NOT consider the sign of gmp_rep
        uint64_t res = mpz_tdiv_ui( (mpz_srcptr)&gmp_rep, l);
#  ifdef __GIVARO_DEBUG
        Integer toto = res;
        if (isneg) toto = -(int64_t)res ;

        // std::cout << toto << ',' << l << ',' << ',' << *this << std::endl;
        assert((toto<l) && (-toto<l) && (toto.priv_sign()*(*this).priv_sign()>=0)) ;
#  endif
        if (!res)
            return (int64_t)res ;
        if (isneg) return (-(int64_t)res) ;


        return (int64_t)res ;
#endif
    }

    int64_t Integer::operator % (const int64_t l) const
    {
        int64_t res ;
        if (l>0) {
            res = static_cast<int64_t>(this->operator%( static_cast<uint64_t>(l) ) );
        }
        else {
            res = static_cast<int64_t>(this->operator%( static_cast<uint64_t>( -l ) ) );
        }

        // std::cout << res << ',' << l << ',' << *this << std::endl;
        assert((res<GIVABS(l)) && (res> -GIVABS(l)) && (((res>0)?1:((res==0)?0:-1))*(*this).priv_sign()>=0)) ;
        return res;
    }

    double Integer::operator % (const double l) const
    {
        double res ;
        if (l>0)
            res =  static_cast<double>(this->operator%( static_cast<uint64_t>(l) ) );
        else{
            res =  static_cast<double>(this->operator%( static_cast<uint64_t>(-l) ) );
        }
        assert((res<GIVABS(l)) && (res> -GIVABS(l)) && (((res>0)?1:((res==0)?0:-1))*(*this).priv_sign()>=0)) ;
        return res;
    }

    // -- operator %
    Integer operator % (const int32_t l, const Integer& n)
    {
        return Integer(l) % n;
    }
    Integer operator % (const int64_t l, const Integer& n)
    {
        return Integer(l) % n;
    }

    Integer operator % (const uint32_t l, const Integer& n)
    {
        return Integer(l) % n;
    }
    Integer operator % (const uint64_t l, const Integer& n)
    {
        return Integer(l) % n;
    }


}


#endif // __GIVARO_gmpxx_gmpxx_int_mod_C

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
