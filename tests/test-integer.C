// Copyright(c)'1994-2011 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// written by BB
// see the COPYRIGHT file for more details.

/*! @file tests/test-integer.C
 * @ingroup tests
 * tests integer.h fucntions not tested elsewhere.
 * @test tests integer.h fucntions not tested elsewhere.
 */

#include <stdlib.h>
#include "gmp++/gmp++.h"

using namespace Givaro;
#define _GIV_REPET 10

int test_axpy_unit(const Integer & a, const Integer & b, const Integer & x)
{
    Integer y,z ;
    /*  axpy */
    z = a*x + b ;
    Integer::axpy(y,a,x,b);
    if (y != z)
        return -1 ;

    y = b ;
    Integer::axpyin(y,a,x);
    if (y != z)
        return -2 ;

    y = b ;
    Integer::axpy(y,a,x,y);
    if (y != z)
        return -3 ;

    /* axmy */
    z = a*x - b ;
    Integer::axmy(y,a,x,b);
    if (y != z)
        return -4 ;

    y = b ;
    Integer::axmyin(y,a,x);
    if (y != z)
        return -5 ;

    y = b ;
    Integer::axmy(y,a,x,y);
    if (y != z)
        return -6 ;

    /* maxpy */
    z =  b - a * x ;
    Integer::maxpy(y,a,x,b);
    if (y != z)
        return -7 ;

    y = b ;
    Integer::maxpyin(y,a,x);
    if (y != z)
        return -8 ;

    y = b ;
    Integer::maxpy(y,a,x,y);
    if (y != z)
        return -9 ;

    return 0 ;

}

int test_divmod()
{
    Integer n, d, q, r;
    uint64_t ud, ur;
    int64_t sd, sr;

    int repet = _GIV_REPET ;
    while (--repet) {
        n = Integer::random<false>();
        d = Givaro::abs(Integer::random<false>());

        // Integer
        q = 0;
        r = 0;
        Integer::divmod(q, r, n, d);
        if (n != q * d + r || r < 0 || r >= d) {
            return -1;
        }

        // uint64_t
        q = 0;
        ur = 0;
        ud = d;
        Integer::divmod(q, ur, n, ud);
        if (n != q * ud + ur || ur >= ud) {
            return -2;
        }

        // int64_t
        q = 0;
        sr = 0;
        sd = d;
        Integer::divmod(q, sr, n, sd);
        if (n != q * sd + sr || sr < 0 || sr >= sd) {
            return -3;
        }
    }

    return 0;
}

int test_axpy()
{
    Integer a,x,b ;
    int res = 0 ;
    int repet = _GIV_REPET ;
    while (--repet) {
        /* axpy */
        a = Integer::random<false>();
        b = Integer::random<false>();
        x = Integer::random<false>();
        res = test_axpy_unit(a,b,x);
        if (res)
            return res ;

        a = 0 ;
        res = test_axpy_unit(a,b,x);
        if (res)
            return res ;

        a = Integer::random<false>();
        b = 0 ;
        res = test_axpy_unit(a,b,x);
        if (res)
            return res ;

        b = Integer::random<false>();
        x = 0 ;
        res = test_axpy_unit(a,b,x);
        if (res)
            return res ;
    }

    return res ;
}

int test_mul_unit(const Integer & a, const Integer & b)
{
    Integer y,z ;
    /*  axpy */
    z = a*b ;
    Integer::mul(y,a,b);
    if (y != z)
        return -1 ;

    y = b ;
    Integer::mulin(y,a);
    if (y != z)
        return -2 ;

    y = a ;
    Integer::mul(y,y,b);
    if (y != z)
        return -3 ;

    y = b ;
    Integer::mul(y,a,y);
    if (y != z)
        return -3 ;

    z = b * b ;
    y = b ;
    Integer::mulin(y,y);
    if (y != z)
        return -5 ;

    return 0 ;

}

int test_mul()
{
    Integer a,b ;
    int res = 0 ;
    int repet = _GIV_REPET ;
    while (--repet) {
        /* axpy */
        a = Integer::random<false>();
        b = Integer::random<false>();
        res = test_mul_unit(a,b);
        if (res)
            return res ;

        a = 0 ;
        res = test_mul_unit(a,b);
        if (res)
            return res ;

        a = Integer::random<false>();
        b = 0 ;
        res = test_mul_unit(a,b);
        if (res)
            return res ;

    }

    return res ;
}

#include <typeinfo>
#include <iostream>


template<typename UnsignedBaseType>
int test_cast_unit(const UnsignedBaseType& t, const Integer& a, const Integer& b) {
#ifdef __GIVARO_DEBUG
    std::cerr << "TYPE: " << typeid(t).name() << std::endl;
#endif

    Integer it(t), iu=t, iv=Integer(t);

    if ( (UnsignedBaseType)it != t) return -1;
    if ( (UnsignedBaseType)iu != t) return -2;
    if ( (UnsignedBaseType)iv != t) return -3;

    UnsignedBaseType at( (UnsignedBaseType)a ), bt( (UnsignedBaseType) b );

    if ( (UnsignedBaseType) (at * bt) != (UnsignedBaseType)(a*b) ) {
#ifdef __GIVARO_DEBUG
        std::cerr << "a: " << a << std::endl;
        std::cerr << "b: " << b << std::endl;
        std::cerr << "a*b: " << (a*b) << std::endl;
        std::cerr << "(a*b)t: " << ( (UnsignedBaseType)(a*b) ) << std::endl;
        std::cerr << "at: " << at << std::endl;
        std::cerr << "bt: " << bt << std::endl;
        std::cerr << "at*bt: " << ( at * bt ) << std::endl;
        std::cerr << "(at*bt)t: " << (UnsignedBaseType)( at * bt ) << std::endl;
#endif
        return -10;
    }

    return 0;
}

template<>
int test_cast_unit(const bool& t, const Integer& a, const Integer& b) {
#ifdef __GIVARO_DEBUG
    std::cerr << "TYPE: " << typeid(t).name() << std::endl;
#endif

    Integer it(t), iu=t, iv=Integer(t);

    if ( (bool)it != t) return -1;
    if ( (bool)iu != t) return -2;
    if ( (bool)iv != t) return -3;

    bool at( (bool)a ), bt( (bool) b );

    if ( (at & bt) != (bool)(a*b) ) {
#ifdef __GIVARO_DEBUG
        std::cerr << "a: " << a << std::endl;
        std::cerr << "b: " << b << std::endl;
        std::cerr << "a*b: " << (a*b) << std::endl;
        std::cerr << "(a*b)t: " << ( (bool)(a*b) ) << std::endl;
        std::cerr << "at: " << at << std::endl;
        std::cerr << "bt: " << bt << std::endl;
        std::cerr << "at&bt: " << ( at & bt ) << std::endl;
#endif
        return -10;
    }

    return 0;
}

int test_cast_unit(const Integer& a, const Integer& b) {
    int res = 0;
    res = test_cast_unit( (bool) Integer::random<false>(), a, b);
    if (res) return res ;
    res = test_cast_unit( (uint16_t) Integer::random<false>(), a, b);
    if (res) return res ;
    res = test_cast_unit( (unsigned char) Integer::random<false>(), a, b);
    if (res) return res ;
    res = test_cast_unit( (uint32_t) Integer::random<false>(), a, b);
    if (res) return res ;
    res = test_cast_unit( (uint64_t) Integer::random<false>(), a, b);
    if (res) return res ;

    return 0;
}


int test_cast() {
    Integer a,b ;
    int res = 0 ;
    int repet = _GIV_REPET ;
    while (--repet) {
        /* axpy */
        a = Integer::random<false>();
        b = Integer::random<false>();
        res = test_cast_unit(a,b);
        if (res) return res ;

        a = 0 ;
        res = test_cast_unit(a,b);
        if (res) return res ;

        a = Integer::random<false>();
        b = 0 ;
        res = test_cast_unit(a,b);
        if (res) return res ;

    }

    return res ;
}

/* See PR#86 */
int test_subin() {
    Integer a,b;

    a = 0;
    a--;
    if (a + 1 != 0) {
#ifdef __GIVARO_DEBUG
        std::cerr << "ERROR1 a: " << a << ", l: " << 1 << std::endl;
#endif
        return -1;
    }

    a = 0;
    a-= INT64_MIN;
    if (a + INT64_MIN != 0) {
#ifdef __GIVARO_DEBUG
        std::cerr << "ERROR2 a: " << a << ", l: " << INT64_MIN << std::endl;
#endif
        return -2;
    }

    a = 0;
    Integer::subin(a,INT64_MIN);
    if (a + INT64_MIN != 0) {
#ifdef __GIVARO_DEBUG
        std::cerr << "ERROR3 a: " << a << ", l: " << INT64_MIN << std::endl;
#endif
        return -3;
    }

    a = 0;
    b = a - INT64_MIN;
    if (b + INT64_MIN != 0) {
#ifdef __GIVARO_DEBUG
        std::cerr << "ERROR4 b: " << b << ", l: " << INT64_MIN << std::endl;
#endif
        return -4;
    }

    return 0;
}


//! @todo test gcd...

int main (int argc, char ** argv)
{
    uint64_t seed = (uint64_t)(argc>1?(uint64_t)atoi(argv[1]):(uint64_t)BaseTimer::seed ());
#ifdef __GIVARO_DEBUG
    std::cerr << "Seed: " << seed << std::endl;
#endif
    Integer::seeding (seed);

    int res = 0 ;
    res = test_axpy();
    if (res)
        return res ;

    res = test_mul();
    if (res)
        return res ;

    res = test_divmod();
    if (res)
        return res ;

    res = test_cast();
    if (res)
        return res ;

    res = test_subin();
    if (res)
        return res ;

    return res ;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
