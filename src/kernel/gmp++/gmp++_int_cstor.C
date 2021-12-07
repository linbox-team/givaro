// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int_cstor.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: M. Samama, T. Gautier
// $Id: gmp++_int_cstor.C,v 1.4 2009-09-17 14:28:22 jgdumas Exp $
// ==========================================================================
/** @file gmp++/gmp++_int_cstor.C
 * cstoring stuff.
 */

#ifndef __GMPplusplus_CSTOR_C__
#define __GMPplusplus_CSTOR_C__
#include <iostream>
#ifndef __GIVARO_INLINE_ALL
#include "gmp++/gmp++.h"
#endif


namespace Givaro {
    //-----------------------------~Integer()
    Integer::~Integer()
    {
        mpz_clear((mpz_ptr)&gmp_rep) ;
    }

    //-------------------------------Integer(const Integer &n)
    Integer::Integer(const Integer &n)
    {
        mpz_init_set ( (mpz_ptr)&gmp_rep, (mpz_srcptr)&(n.gmp_rep)) ;
    }

    //-----------------------------Integer(int32_t n)
    Integer::Integer(int32_t n)
    {
        mpz_init_set_si((mpz_ptr)&gmp_rep, n) ;
    }

    //-----------------------------Integer(uint32_t n)
    Integer::Integer(unsigned char n)
    {
        mpz_init_set_ui((mpz_ptr)&gmp_rep, n) ;
    }

    //-----------------------------Integer(uint32_t n)
    Integer::Integer(uint32_t n)
    {
        mpz_init_set_ui((mpz_ptr)&gmp_rep, n) ;
    }

    //-----------------------------Integer(int64_t n)
    Integer::Integer(int64_t n)
    {
#if __GIVARO_SIZEOF_LONG < 8
        // 64 bits is less than 20 digits
        char * tmp = new char[23];
        sprintf(tmp,"%lld",n);
        mpz_init_set_str((mpz_ptr)&gmp_rep, tmp, 10) ;
        delete [] tmp;
#else
        mpz_init_set_si((mpz_ptr)&gmp_rep, n) ;
#endif
    }

    //-----------------------------Integer(uint64_t n)
    Integer::Integer(uint64_t n)
    {
#if __GIVARO_SIZEOF_LONG < 8
        // 64 bits is less than 20 digits
        char * tmp = new char[23];
        sprintf(tmp,"%llu",n);
        mpz_init_set_str((mpz_ptr)&gmp_rep, tmp, 10) ;
        delete [] tmp;
#else
        mpz_init_set_ui((mpz_ptr)&gmp_rep, n) ;
#endif
    }

    //-----------------------------Integer(double)
    Integer::Integer(double d)
    {
        mpz_init_set_d((mpz_ptr)&gmp_rep, d) ;
    }

    // -- Integer(const char *s)
    Integer::Integer(const char *s)
    {
        mpz_init_set_str((mpz_ptr)&gmp_rep, s, 10);
    }

    //------------------------------------------operator = (const Integer &n)
    Integer& Integer::logcpy(const Integer &n)
    {
        if (this == &n) return *this;
        mpz_set ( (mpz_ptr)&gmp_rep, (mpz_srcptr)&(n.gmp_rep)) ;
        return *this;
    }

    // same as logcopy
    Integer& Integer::operator = (const Integer &n)
    {
        return logcpy(n) ;
    }


    Integer& Integer::copy(const Integer &n)
    {
        if (this == &n) return *this;
        mpz_set ( (mpz_ptr)&gmp_rep, (mpz_srcptr)&(n.gmp_rep)) ;
        return *this ;
    }

    namespace Protected {
        void importWords(Integer& x, size_t count, int32_t order, int32_t size,
                         int32_t endian, size_t nails, const void* op)
        {
            mpz_import( (mpz_ptr)&(x.gmp_rep), count, order, size, endian, nails, op);
        }
    }

    Integer::Integer(const vect_t & v)
    {
        size_t s = v.size();
        if (s) {
            mpz_init_set_ui((mpz_ptr)&gmp_rep, v[0]);
            Integer base(256), prod, tmp;
            prod = base = pow(base, (uint64_t)sizeof(mp_limb_t) );

            std::vector<mp_limb_t>::const_iterator vi = v.begin();
            for(++vi;vi != v.end();++vi) {
                mpz_mul_ui( (mpz_ptr)&tmp.gmp_rep, (mpz_ptr)&prod.gmp_rep, *vi);
                *this += tmp;
                prod *= base;
            }
        } else
            mpz_init( (mpz_ptr)&gmp_rep );

    }


}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
