// ========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: A. Breust
// Time-stamp: <23 Sep 16 14:03:02 Jean-Guillaume.Dumas@imag.fr>
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
};
template<> struct IntType<true> {
    typedef int32_t type;
};

template<typename T>
struct make_signed_int<T,
    typename std::enable_if<std::is_floating_point<T>::value>::type> {
    typedef typename IntType<std::is_same<T,float>::value>::type type;
};

namespace Givaro
{

    template<typename Storage_t>
    inline Storage_t& invext(Storage_t& x, Storage_t& d, const Storage_t a, const Storage_t b)
    {
        using Compute_t = typename make_signed_int<Storage_t>::type;

        Compute_t u0(0), r0 = static_cast<Compute_t>(b);
        Compute_t u1(1), r1 = static_cast<Compute_t>(a);
        while (r1 != static_cast<Compute_t>(0))
        {
            Compute_t q(static_cast<Compute_t>(r0 / r1));
            Compute_t t(u1);
            u1 = static_cast<Compute_t>(u0 - q * u1);
            u0 = t;

            t = r1;
            r1 = static_cast<Compute_t>(r0 - q * r1);
            r0 = t;
//         std::cerr << u0 << '*' << a << "+ .*" << b << '=' << r0 << std::endl;
//         std::cerr << u1 << '*' << a << "+ .*" << b << '=' << r1 << std::endl;
        }
        d = static_cast<Storage_t>(r0);
//         std::cerr << u0 << '*' << a << "+ .*" << b << '=' << r0 << std::endl;
        return x = ( u0<0 ? static_cast<Storage_t>(u0+b) : static_cast<Storage_t>(u0) );
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
        invext(u,d,a,b);
        v = static_cast<Storage_t>((d-u*a)/b);
        return d;
    }

    template<>
    inline float invext(const float a, const float b)
    {
        float u1(1.f), v1(0.f);
	float u3(a), v3(b);

        while (v3 != 0.f)
        {
	    float q = std::floor(u3 / v3);
            float t;
	    
	    t = v1;
            v1 = u1 - q * v1;
            u1 = t;

            t = v3;
            v3 = u3 - q * v3;
            u3 = t;
        }

        return u1;
    }

    template<>
    inline double invext(const double a, const double b)
    {
        double u1(1.0), v1(0.0);
	double u3(a), v3(b);

        while (v3 != 0.0)
        {
	    double q = std::floor(u3 / v3);
            double t;
	    
	    t = v1;
            v1 = u1 - q * v1;
            u1 = t;

            t = v3;
            v3 = u3 - q * v3;
            u3 = t;
        }

        return u1;
    }

}

