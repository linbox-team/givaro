// ========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: A. Breust
// Time-stamp: <28 May 15 09:40:00 Alexis.Breust@imag.fr>
// ========================================================================
// Description:
// Forward declarations for Givaro::Modular and associated functions

#pragma once

namespace Givaro
{
    template<typename Storage_t>
    inline Storage_t& gcdext(Storage_t& d, Storage_t& u, Storage_t& v, const Storage_t a, const Storage_t b)
    {
        using Compute_t = typename std::make_signed<Storage_t>::type;

        Compute_t u1 = 1, u2 = 0, u3 = static_cast<Compute_t>(a);
        Compute_t v1 = 0, v2 = 1, v3 = static_cast<Compute_t>(b);

        while (v3 != 0)
        {
            Compute_t q = static_cast<Compute_t>(u3 / v3);
            Compute_t t1, t2 , t3;

            // Weirdly, all the casts are needed when used with uint8_t,
            // as gcc operators on this type returns an int, not a uint8_t.
            t1 = static_cast<Compute_t>(u1 - q * v1);
            t2 = static_cast<Compute_t>(u2 - q * v2);
            t3 = static_cast<Compute_t>(u3 - q * v3);

            u1 = v1;
            u2 = v2;
            u3 = v3;

            v1 = t1;
            v2 = t2;
            v3 = t3;
        }

        u = static_cast<Storage_t>(u1);
        v = static_cast<Storage_t>(u2);
        return d = static_cast<Storage_t>(u3);
    }

    template<typename Storage_t>
    inline Storage_t& invext(Storage_t& x, const Storage_t a, const Storage_t b)
    {
        using Compute_t = typename std::make_signed<Storage_t>::type;

        Compute_t u1(1), u3 = static_cast<Compute_t>(a);
        Compute_t v1(0), v3 = static_cast<Compute_t>(b);

        while (v3 != static_cast<Compute_t>(0))
        {
            Compute_t q = static_cast<Compute_t>(u3 / v3);
            Compute_t t;

            t = v1;
            v1 = static_cast<Compute_t>(u1 - q * v1);
            u1 = t;

            t = v3;
            v3 = static_cast<Compute_t>(u3 - q * v3);
            u3 = t;
        }

        return x = static_cast<Storage_t>(u1);
    }

    template<typename Storage_t>
    inline Storage_t invext(const Storage_t a, const Storage_t b)
    {
        Storage_t r;
        return invext(r, a, b);
    }
}

