// ========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: A. Breust
// Time-stamp: <12 May 21 09:46:39 Jean-Guillaume.Dumas@imag.fr>
// ========================================================================
// Description:
// Forward declarations for Givaro::Modular and associated functions

#pragma once

#include <cmath>
#include <type_traits>

// int64_t for double and int32_t for float
template<bool IsFloat=false> struct IntType {
    typedef int64_t type;
    typedef uint64_t utype;
};
template<> struct IntType<true> {
    typedef int32_t type;
    typedef uint32_t utype;
};

template<typename T, typename Enable = void>
struct make_unsigned_int {
    typedef typename std::make_unsigned<T>::type type;
};
template<typename T, typename Enable = void>
struct make_signed_int {
    typedef typename std::make_signed<T>::type type;
};

template<typename T>
struct make_unsigned_int<T,
typename std::enable_if<std::is_floating_point<T>::value>::type> {
    typedef typename IntType<std::is_same<T,float>::value>::utype type;
};

template<typename T>
struct make_signed_int<T,
typename std::enable_if<std::is_floating_point<T>::value>::type> {
    typedef typename IntType<std::is_same<T,float>::value>::type type;
};

template<typename> struct is_ruint : std::false_type {};
template<size_t K> struct is_ruint<RecInt::ruint<K>> : std::true_type {};

template<typename, typename> struct is_same_ruint : std::false_type {};
template<size_t K> struct is_same_ruint<RecInt::ruint<K>,RecInt::ruint<K>> : std::true_type {};

template<typename, typename> struct is_smaller_ruint : std::false_type {};
template<size_t K> struct is_smaller_ruint<RecInt::ruint<K>,RecInt::ruint<K+1>> : std::true_type {};

template<typename, typename> struct is_smaller_rint : std::false_type {};
template<size_t K> struct is_smaller_rint<RecInt::rint<K>,RecInt::rint<K+1>> : std::true_type {};

template<typename T> struct RecInt_K{ static const size_t value = 0;};
template<size_t K> struct RecInt_K<RecInt::ruint<K>> { static const size_t value = K;};
template<size_t K> struct RecInt_K<RecInt::rint<K>>  { static const size_t value = K;};


template<typename> struct is_rint : std::false_type {};
template<size_t K> struct is_rint<RecInt::rint<K>> : std::true_type {};

template<typename, typename> struct is_same_rint : std::false_type {};
template<size_t K> struct is_same_rint<RecInt::rint<K>,RecInt::rint<K>> : std::true_type {};


namespace Givaro
{

    template<typename Storage_t>
    inline typename std::enable_if<!std::is_floating_point<Storage_t>::value, Storage_t&>::type
    extended_euclid (Storage_t& x, Storage_t& d, const Storage_t a, const Storage_t b)
    {
        Storage_t u0(0), u1(1), r1(a), q, t;
        d = b;

        bool neg = true;

        while (r1 != (Storage_t)0)
        {
            // Invariants:
            // - if !neg:
            //     u0 * a = d (mod b)
            //     u1 * a = -r1 (mod b)
            // - if neg:
            //     u0 * a = -d (mod b)
            //     u1 * a = r1 (mod b)

            q = d / r1;
            t = u1;
            u1 = q * u1 + u0;
            u0 = t;

            t = r1;
            r1 = d - q * r1;
            d = t;

            neg = !neg;
        }

        return x = (neg && u0 > 0) ? b - u0 : u0;
    }

    template<typename Storage_t>
    inline typename std::enable_if<std::is_floating_point<Storage_t>::value, Storage_t&>::type
    extended_euclid(Storage_t& x, Storage_t& d, const Storage_t a, const Storage_t b)
    {
        Storage_t u1(1), v1(0);
        Storage_t u3(a), v3(b);

        while (v3 != 0)
        {
            Storage_t q = std::floor(u3 / v3);
            Storage_t t;

            t = v1;
            v1 = u1 - q * v1;
            u1 = t;

            t = v3;
            v3 = u3 - q * v3;
            u3 = t;
        }

        d = u3;

        return x = u1;
    }

    template<typename Storage_t>
    inline Storage_t& invext (Storage_t& x, Storage_t& d, const Storage_t a, const Storage_t b){

        extended_euclid (x,d,a,b);
#ifdef __GIVARO_DEBUG
        if ( d > (Storage_t)1 ) {
            throw GivMathDivZero("*** Error: division by zero, in operator invext<floating_point> in modular-general") ;
        }
#endif
        return x;
    }

    template<typename Storage_t>
    inline Storage_t& invext(Storage_t& x, const Storage_t a, const Storage_t b)
    {
        Storage_t g; return invext(x,g,a,b);
    }

    template<typename Storage_t>
    inline Storage_t invext(const Storage_t a, const Storage_t b)
    {
        Storage_t r; return invext(r, a, b);
    }

    template<typename Storage_t>
    inline Storage_t& gcdext(Storage_t& d, Storage_t& u, Storage_t& v, const Storage_t a, const Storage_t b)
    {
        extended_euclid (u,d,a,b);
        v = (d-u*a)/b;
        return d;
    }


}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
