// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int_misc.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: M. Samama, T. Gautier
// $Id: gmp++_int_misc.C,v 1.16 2010-12-16 16:54:38 jgdumas Exp $
// ==========================================================================
// Description:
/** @file gmp++/gmp++_int_misc.C
 * miscing stuff.
 */

#ifndef __GIVARO_gmpxx_gmpxx_int_misc_C
#define __GIVARO_gmpxx_gmpxx_int_misc_C

#include <iostream>
#include <cmath>
#ifndef __GIVARO_INLINE_ALL
#include "gmp++/gmp++.h"
#endif
#ifndef __GIVARO_GMP_NO_CXX
#include <sstream>
#endif

namespace Givaro {
    //-------------------------------------------fact (uint64_t l)
    Integer fact ( uint64_t l)
    {
        Integer Res ;
#if __GIVARO_SIZEOF_LONG < 8
        // factorial too large if l > 2^32 anyway
        uint32_t i(l);
        GIVARO_ASSERT( (uint64_t)i == l, "Factorial too large");
        mpz_fac_ui( (mpz_ptr)&(Res.gmp_rep), i ) ;
#else
        mpz_fac_ui( (mpz_ptr)&(Res.gmp_rep), l ) ;
#endif
        return Res ;
    }

    //-------------------------------------------square root
    Integer& sqrt(Integer& q, const Integer &a)
    {
        mpz_sqrt( (mpz_ptr)&(q.gmp_rep),
                  (mpz_srcptr)&(a.gmp_rep)) ;
        return q;
    }

    Integer& sqrtrem(Integer& q, const Integer &a, Integer& r)
    {
        mpz_sqrtrem( (mpz_ptr)&(q.gmp_rep),
                     (mpz_ptr)&(r.gmp_rep), (mpz_srcptr)&(a.gmp_rep)) ;
        return q;
    }

    Integer sqrt(const Integer &a)
    {
        Integer q;
        return sqrt(q,a);
    }

    Integer sqrtrem(const Integer &a, Integer& r)
    {
        Integer q;
        return sqrtrem(q,a,r);
    }

    bool root(Integer& q, const Integer &a, uint32_t n)
    {
        return (bool)mpz_root ((mpz_ptr)&(q.gmp_rep),
                               (mpz_srcptr)&(a.gmp_rep),
                               n);
    }

    void swap(Integer& a, Integer& b)
    {
        return mpz_swap( (mpz_ptr)&(a.gmp_rep), (mpz_ptr)&(b.gmp_rep));
    }


    // Natural logarithm of a
    // log(2) being close to 0.69314718055994531
    double naturallog(const Integer& a)
    {
        long int exp_;
        double d = mpz_get_d_2exp( &exp_, (mpz_srcptr)&(a.gmp_rep) );
        return (double)exp_*0.69314718055994531+log(d);
    }


    /*! Tests parity of an integer
     * @param a integer
     * @return 1 if odd, 0 if even
     */
    bool isOdd(const Integer &a)
    {
        int32_t o = mpz_tstbit( (mpz_srcptr) &(a.gmp_rep), 0);
        return (o!=0); // or maybe should I write l==1 ^^
    }

    // base p logarithm of a
    int64_t logp(const Integer& a, const Integer& p)
    {
        std::list< Integer > pows;
        Integer puiss = p, sq;
        do {
            pows.push_back( puiss );
        } while ( (puiss *= puiss) <= a );
        puiss = pows.back(); pows.pop_back();
        int64_t res = (1 << pows.size());
        while (! pows.empty() ) {
            if ((sq = puiss * pows.back()) <= a) {
                puiss = sq;
                pows.pop_back();
                res += (1 << pows.size());
            } else
                pows.pop_back();
        }
        return res;
    }

    // approximation of the base 2 logarithm of a
    // 1/log(2) being close to 1.44269504088896341
    double logtwo(const Integer& a)
    {
        long int exp;
        double d = mpz_get_d_2exp( &exp, (mpz_srcptr)&(a.gmp_rep) );
        return (double)exp+log(d)*1.44269504088896341;
    }

    namespace Protected {
        //------------------------------------------GMP isprime
        //     If this function returns 0, OP is definitely not prime.  If it
        //     returns 1, then OP is `probably' prime.  The probability of a
        //     false positive is (1/4)^r.  A reasonable value of r is 25.
        int32_t probab_prime(const Integer &p, int32_t r)
        {
            return mpz_probab_prime_p ((mpz_srcptr)&(p.gmp_rep),r) ;
        }

        Integer& nextprime(Integer& r, const Integer &p)
        {
            mpz_nextprime ((mpz_ptr)&(r.gmp_rep), (mpz_srcptr)&(p.gmp_rep)) ;
            return r;
        }

        // Copied and adapted from mpz/nextprime.c
        Integer& prevprime(Integer& r, const Integer &p)
        {
            if (p < 3) return (r=2u);
            if (isOdd(p))
                mpz_sub_ui ( (mpz_ptr)&(r.gmp_rep), (mpz_srcptr)&(p.gmp_rep), 2u );
            else
                mpz_sub_ui ( (mpz_ptr)&(r.gmp_rep), (mpz_srcptr)&(p.gmp_rep), 1u );

            while( !mpz_probab_prime_p ( (mpz_srcptr)&(r.gmp_rep), _GIVARO_ISPRIMETESTS_ ) )
            {

                mpz_sub_ui ( (mpz_ptr)&(r.gmp_rep), (mpz_srcptr)&(r.gmp_rep), 2u );

            }

            return r;
        }
    } // namespace Protected


    // ==========================================================================
    // Computes and returns the Kronecker/Jacobi and Legendre symbols (u/v)
    // of the integers u and v.
    // The algorithm used is Gmp's.
    int32_t kronecker(const Integer& u, const Integer& v)
    {
        return mpz_kronecker ((mpz_srcptr)&(u.gmp_rep),(mpz_srcptr)&(v.gmp_rep)) ;
    }

    int32_t jacobi(const Integer& u, const Integer& v)
    {
        return mpz_jacobi ((mpz_srcptr)&(u.gmp_rep),(mpz_srcptr)&(v.gmp_rep)) ;
    }

    int32_t legendre(const Integer& u, const Integer& v)
    {
        return mpz_legendre ((mpz_srcptr)&(u.gmp_rep),(mpz_srcptr)&(v.gmp_rep)) ;
    }



    //--------------------------------------------Integer::operator <<   // shift left
    Integer Integer::operator << (int32_t l) const
    {
        return this->operator<<( (uint64_t)l );
    }
    Integer Integer::operator << (uint32_t l) const
    {
        return this->operator<<( (uint64_t)l );
    }
    Integer Integer::operator << (int64_t l) const
    {
        return this->operator<<( (uint64_t)l );
    }

    Integer Integer::operator << (uint64_t l) const
    {
        Integer tmp;
        mpz_mul_2exp((mpz_ptr)&(tmp.gmp_rep), (mpz_srcptr)&(gmp_rep), l );
        return tmp;
    }


    //--------------------------------------------Integer::operator >>   // shift right
    Integer Integer::operator >> (int32_t l) const
    {
        return this->operator>>( (uint64_t)l );
    }

    Integer Integer::operator >> (int64_t l) const
    {
        return this->operator>>( (uint64_t)l );
    }

    Integer Integer::operator >> (uint32_t l) const
    {
        return this->operator>>( (uint64_t)l );
    }

    Integer Integer::operator >> (uint64_t l) const
    {
        Integer tmp;
        mpz_tdiv_q_2exp( (mpz_ptr)&(tmp.gmp_rep), (mpz_srcptr)&(gmp_rep), l );
        return tmp;
    }

    //--------------------------------------------Integer::operator <<=   // shift left
    Integer& Integer::operator <<= (int32_t l)
    {
        return this->operator<<= ( (uint64_t)l );
    }
    Integer& Integer::operator <<=  (uint32_t l)
    {
        return this->operator<<= ( (uint64_t)l );
    }
    Integer& Integer::operator <<= (int64_t l)
    {
        return this->operator<<= ( (uint64_t)l );
    }

    Integer& Integer::operator <<= (uint64_t l)
    {
        mpz_mul_2exp((mpz_ptr)&(gmp_rep), (mpz_srcptr)&(gmp_rep), l );
        return *this;
    }


    //--------------------------------------------Integer::operator >>=   // shift right
    Integer& Integer::operator >>= (int32_t l)
    {
        return this->operator>>= ( (uint64_t)l );
    }
    Integer& Integer::operator >>= (int64_t l)
    {
        return this->operator>>= ( (uint64_t)l );
    }
    Integer& Integer::operator >>= (uint32_t l)
    {
        return this->operator>>= ( (uint64_t)l );
    }

    Integer& Integer::operator >>= (uint64_t l)
    {
        mpz_tdiv_q_2exp( (mpz_ptr)&(gmp_rep), (mpz_srcptr)&(gmp_rep), l );
        return *this;
    }

    //------------------------------------------- Bit logic
    Integer Integer::operator^ (const Integer& a) const
    {   // XOR
        Integer res(*this);
        return res ^= a;
    }
    Integer Integer::operator| (const Integer& a) const
    {   // OR
        Integer res(*this);
        return res |= a;
    }
    Integer Integer::operator& (const Integer& a) const
    {   // AND
        Integer res(*this);
        return res &= a;
    }
    Integer Integer::operator^ (const uint64_t & a) const
    {   // XOR
        Integer res(*this);
        return res ^= a;
    }
    Integer Integer::operator| (const uint64_t & a) const
    {   // OR
        Integer res(*this);
        return res |= a;
    }
    uint64_t Integer::operator& (const uint64_t & a) const
    {   // AND
        return mpz_get_ui((mpz_srcptr)&(gmp_rep)) & a;
    }
    Integer Integer::operator^ (const uint32_t& a) const
    {   // XOR
        Integer res(*this);
        return res ^= a;
    }
    Integer Integer::operator| (const uint32_t& a) const
    {   // OR
        Integer res(*this);
        return res |= a;
    }
    uint32_t Integer::operator& (const uint32_t& a) const
    {   // AND
        return (uint32_t) (mpz_get_ui((mpz_srcptr)&(gmp_rep)) & (uint64_t)a );
    }
    Integer Integer::operator~ () const
    {   // 1 complement
        Integer res;
        mpz_com( (mpz_ptr)&(res.gmp_rep), (mpz_srcptr)&(gmp_rep));
        return res;
    }
    Integer& Integer::operator^= (const Integer& a)
    {   // XOR
        mpz_xor( (mpz_ptr)&(gmp_rep), (mpz_ptr)&(gmp_rep), (mpz_srcptr)&(a.gmp_rep));
        return *this;
    }
    Integer& Integer::operator|= (const Integer& a)
    {   // OR
        mpz_ior( (mpz_ptr)&(gmp_rep), (mpz_ptr)&(gmp_rep), (mpz_srcptr)&(a.gmp_rep));
        return *this;
    }
    Integer& Integer::operator&= (const Integer& a)
    {   // AND
        mpz_and( (mpz_ptr)&(gmp_rep), (mpz_ptr)&(gmp_rep), (mpz_srcptr)&(a.gmp_rep));
        return *this;
    }

    Integer& Integer::operator^= (const uint64_t & a)
    {   // XOR
        Integer au(a);
        mpz_xor( (mpz_ptr)&(gmp_rep), (mpz_ptr)&(gmp_rep), (mpz_srcptr)&(au.gmp_rep));
        return *this;
    }
    Integer& Integer::operator|= (const uint64_t & a)
    {   // OR
        Integer au(a);
        mpz_ior( (mpz_ptr)&(gmp_rep), (mpz_ptr)&(gmp_rep), (mpz_srcptr)&(au.gmp_rep));
        return *this;
    }
    Integer& Integer::operator&= (const uint64_t & a)
    {   // AND
        Integer au(a);
        mpz_and( (mpz_ptr)&(gmp_rep), (mpz_ptr)&(gmp_rep), (mpz_srcptr)&(au.gmp_rep));
        return *this;
    }

    Integer& Integer::operator^= (const uint32_t& a)
    {   // XOR
        Integer au(a);
        mpz_xor( (mpz_ptr)&(gmp_rep), (mpz_ptr)&(gmp_rep), (mpz_srcptr)&(au.gmp_rep));
        return *this;
    }
    Integer& Integer::operator|= (const uint32_t& a)
    {   // OR
        Integer au(a);
        mpz_ior( (mpz_ptr)&(gmp_rep), (mpz_ptr)&(gmp_rep), (mpz_srcptr)&(au.gmp_rep));
        return *this;
    }
    Integer& Integer::operator&= (const uint32_t& a)
    {   // AND
        Integer au(a);
        mpz_and( (mpz_ptr)&(gmp_rep), (mpz_ptr)&(gmp_rep), (mpz_srcptr)&(au.gmp_rep));
        return *this;
    }


    //------------------------------------------- convert method
    //------------------------------------------- casting method
    Integer::operator int32_t() const
    {
        return int32_t (mpz_get_si ( (mpz_srcptr)&gmp_rep));
    }
    Integer::operator uint32_t() const
    {
        return (uint32_t) mpz_get_ui ( (mpz_srcptr)&gmp_rep);
    }
    Integer::operator int64_t() const
    {
#if __GIVARO_SIZEOF_LONG < 8
        uint64_t u64 = this->operator uint64_t();
        /* The following lines are adapted from the file mpz/get_si.c of GMP
         * source code.
         */
        auto sgn = this->sign();
        if (sgn > 0)
            return u64 & INT64_MAX;
        else if (sgn < 0)
            return -1 - (int64_t) ((u64 -1) & INT64_MAX);
        else
            return 0;
#else
        return mpz_get_si ( (mpz_srcptr)&gmp_rep);
#endif
    }
    Integer::operator uint64_t() const
    {
#if __GIVARO_SIZEOF_LONG < 8
        Integer absThis = abs(*this);
        uint64_t r = static_cast<uint64_t>(absThis.operator uint32_t());
        absThis >>= 32;
        return r |= static_cast<uint64_t>(absThis.operator uint32_t()) << 32;
#else
        return mpz_get_ui ( (mpz_srcptr)&gmp_rep);
#endif
    }
    Integer::operator double() const
    {
        return mpz_get_d ( (mpz_srcptr)&gmp_rep);
    }
    Integer::operator float() const
    {
        return (float)mpz_get_d ( (mpz_srcptr)&gmp_rep);
    }

    Integer::operator std::string () const
    {
        std::ostringstream o ;
        print(o);
        return o.str();
    }

    Integer::operator Integer::vect_t () const
    {
        size_t s = mpz_size( (mpz_srcptr)&(gmp_rep) );
        std::vector<mp_limb_t> v(s);
        std::vector<mp_limb_t>::iterator vi = v.begin();
        for(mp_size_t i = 0;vi != v.end();++vi, ++i)
            *vi = mpz_getlimbn( (mpz_srcptr)& (gmp_rep) ,i);
        return v;
    }

    uint64_t length(const Integer& a)
    {
        //! @bug JGD 23.04.2012: shouldn't it be "mp_limb_t" instead of "uint64_t"?
        return mpz_size( (mpz_srcptr)&(a.gmp_rep) ) * sizeof(uint64_t);
    }

    Integer abs(const Integer &n)
    {
        if (sign(n) >= 0)
            return n;
        return -n;
    }

    size_t Integer::size() const
    {
        return  mpz_size( (mpz_srcptr)&gmp_rep ) ;
    }

    size_t Integer::size_in_base(int32_t BASE) const
    {
        return  mpz_sizeinbase ((mpz_srcptr)&gmp_rep, BASE);
    }

    size_t Integer::bitsize() const
    {
        return  mpz_sizeinbase ((mpz_srcptr)&gmp_rep, 2);
    }

    uint64_t Integer::operator[](size_t i) const
    {
        if ( mpz_size( (mpz_srcptr)&gmp_rep ) > i)
            return mpz_getlimbn( (mpz_srcptr)&gmp_rep, i);
        else
            return 0;
    }

}

#endif // __GIVARO_gmpxx_gmpxx_int_misc_C

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
