// ========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: A. Breust
// Time-stamp: <27 Sep 16 18:53:18 Jean-Guillaume.Dumas@imag.fr>
// ========================================================================
// Description:
// Forward declarations for Givaro::Modular and associated functions

#pragma once

#include <cmath>
#include <type_traits>

template<typename T, typename Enable = void>
struct make_signed_int {
    typedef typename std::make_signed<T>::type type;
};

// int64_t for double and int32_t for float
template<bool IsFloat=false> struct IntType {
    typedef int64_t type;
    typedef uint64_t utype;
};
template<> struct IntType<true> {
    typedef int32_t type;
    typedef uint32_t utype;
};

template<typename T>
struct make_signed_int<T,
    typename std::enable_if<std::is_floating_point<T>::value>::type> {
    typedef typename IntType<std::is_same<T,float>::value>::type type;
};

template<typename T, typename Enable = void>
struct make_unsigned_int {
    typedef typename std::make_unsigned<T>::type type;
};

template<typename T>
struct make_unsigned_int<T,
    typename std::enable_if<std::is_floating_point<T>::value>::type> {
    typedef typename IntType<std::is_same<T,float>::value>::utype type;
};

template<typename> struct is_RecInt : std::false_type {};
template<size_t K> struct is_RecInt<RecInt::ruint<K>> : std::true_type {};
template<size_t K> struct is_RecInt<RecInt::rint<K>> : std::true_type {};

template<typename, typename> struct is_same_RecInt : std::false_type {};
template<size_t K> struct is_same_RecInt<RecInt::ruint<K>,RecInt::ruint<K>> : std::true_type {};
template<size_t K> struct is_same_RecInt<RecInt::rint<K>,RecInt::rint<K>> : std::true_type {};

template<typename, typename> struct is_smaller_RecInt : std::false_type {};
template<size_t K> struct is_smaller_RecInt<RecInt::ruint<K>,RecInt::ruint<K+1>> : std::true_type {};

template<typename T> struct RecInt_K;
template<size_t K> struct RecInt_K<RecInt::ruint<K>> { static const size_t value = K;};



namespace Givaro
{

    template<typename Storage_t, typename Compute_t>
    //inline Storage_t& 
    inline typename std::enable_if<!std::is_floating_point<Storage_t>::value, Storage_t&>::type
    invext(Storage_t& x, Storage_t& d, const Compute_t a, const Compute_t b)
    {
        //using Compute_t = typename make_signed_int<Storage_t>::type;

        Compute_t u0(0), r0 = b; //Caster<Compute_t>(b);
        Compute_t u1(1), r1 = a; //Caster<Compute_t>(a);
        while (r1 != (Compute_t)0)
        {
            Compute_t q(r0 / r1);
            Compute_t t(u1);
            u1 = u0 - q * u1;
            u0 = t;

            t = r1;
            r1 = r0 - q * r1;
            r0 = t;
//         std::cerr << u0 << '*' << a << "+ .*" << b << '=' << r0 << std::endl;
//         std::cerr << u1 << '*' << a << "+ .*" << b << '=' << r1 << std::endl;
        }
        d = Caster<Storage_t>(r0);
//         std::cerr << u0 << '*' << a << "+ .*" << b << '=' << r0 << std::endl;
        return x = ( u0<0 ? Caster<Storage_t>(u0+Caster<Compute_t>(b)) : Caster<Storage_t>(u0) );
    }

    template<typename Storage_t, typename Compute_t>
    inline typename std::enable_if<std::is_floating_point<Storage_t>::value, Storage_t&>::type
    invext(Storage_t& x, Storage_t& d, const Compute_t a, const Compute_t b)
    {
        Compute_t u1(1), v1(0);
        Compute_t u3(a), v3(b);

        while (v3 != (Compute_t)0)
        {
            Compute_t q = std::floor(u3 / v3);
            Compute_t t;
            
            t = v1;
            v1 = u1 - q * v1;
            u1 = t;

            t = v3;
            v3 = u3 - q * v3;
            u3 = t;
        }

        d = Caster<Storage_t>(u3); 
        return x = Caster<Storage_t>(u1);
    }

    //template<>
    //inline float invext(const float a, const float b)
    //{
    //    float u1(1.f), v1(0.f);
    //    float u3(a), v3(b);

    //    while (v3 != 0.f)
    //    {
    //        float q = std::floor(u3 / v3);
    //        float t;
    //        
    //        t = v1;
    //        v1 = u1 - q * v1;
    //        u1 = t;

    //        t = v3;
    //        v3 = u3 - q * v3;
    //        u3 = t;
    //    }

    //    return u1;
    //}

    //template<>
    //inline double invext(const double a, const double b)
    //{
    //    double u1(1.0), v1(0.0);
    //    double u3(a), v3(b);

    //    while (v3 != 0.0)
    //    {
    //        double q = std::floor(u3 / v3);
    //        double t;
    //        
    //        t = v1;
    //        v1 = u1 - q * v1;
    //        u1 = t;

    //        t = v3;
    //        v3 = u3 - q * v3;
    //        u3 = t;
    //    }

    //    return u1;
    //}

    template<typename Storage_t, typename Compute_t>
    inline Storage_t& invext(Storage_t& x, const Compute_t a, const Compute_t b)
    {
        Storage_t g; return invext(x,g,a,b);
    }

    template<typename Storage_t, typename Compute_t>
    inline Storage_t invext(const Compute_t a, const Compute_t b)
    {
        Storage_t r; return invext(r, a, b);
    }

    template<typename Storage_t, typename Compute_t>
    inline Storage_t& gcdext(Storage_t& d, Storage_t& u, Storage_t& v, const Compute_t a, const Compute_t b)
    {
        invext(u,d,a,b);
        v = static_cast<Storage_t>((d-u*a)/b);
        return d;
    }


}

