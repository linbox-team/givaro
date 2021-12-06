// ============================================================= //
// Copyright(c)'2011-2016 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors : Jean-Guillaume Dumas (Jean-Guillaume.Dumas@imag.fr)
//           Alexis Breust (alexis.breust@gmail.com 2014)
//           Christophe Chabot (christophechabotcc@gmail.com 2011)
// ============================================================= //
// Function:
// Conversion functions between r(u/m)int and mpz_class from GMP
// ============================================================= //


#ifndef RUINT_CONVERT_H
#define RUINT_CONVERT_H

#include <cstddef> // required by gmp versions <= 5.1.3
#include <gmpxx.h>
#include "rumanip.h" /* reset() */

// --------------------------------------------------------------
// ----------------------- DEFINITIONS ---------------------------

namespace RecInt
{
    // Converts a mpz_class to a ruint<K> and vice-versa
    template <size_t K> ruint<K>& mpz_to_ruint(ruint<K>&, const mpz_class&);
    template <size_t K> mpz_class& ruint_to_mpz(mpz_class&, const ruint<K>&);
}


// --------------------------------------------------------------
// --------------------- Convert ruint --------------------------

namespace RecInt
{
    // Convert a GMP integer into a ruint
    template <size_t K>
    inline ruint<K>& mpz_to_ruint(ruint<K>& a, const mpz_class& b) {
        unsigned int i;
        mpz_class c(b);

        reset(a);
        for (i = 0; i < NBLIMB<K>::value; i++) {
#if __GIVARO_SIZEOF_LONG < 8
            limb l = c.get_ui(); c >>= 32;
            l |= (limb(c.get_ui()) << 32);
            set_limb(a, l, i); c >>= 32;
#else
            limb l = c.get_ui();
            set_limb(a, l, i); c >>= 64;
#endif
        }


        return a;
    }

    // Convert a ruint into a GMP integer
    template <size_t K>
    inline mpz_class& ruint_to_mpz(mpz_class& a, const ruint<K>& b) {
        //         a = 0;
        //         for (auto it(b.rbegin()); it != b.rend(); ++it) {
        // #if __GIVARO_SIZEOF_LONG < 8
        // 		    // GMP does not handle uint64_t, need to break it
        //             a <<= 32;
        //             a ^= static_cast<uint32_t>(mp_limb_t((*it) >> 32));
        //             a <<= 32;
        //             a += static_cast<uint32_t>(mp_limb_t(*it));
        // #else
        //             a <<= 64;
        // #if __GIVARO_SIZEOF_LONG == 8
        //             a ^= static_cast<unsigned long>(*it);
        // #else
        //             a ^= static_cast<uint64_t>(*it);
        // #endif
        // #endif
        //         }
        // uses the fact that recint words will be contiguous
        mpz_import(a.get_mpz_t(), NBLIMB<K>::value, -1, sizeof(limb), 0, 0, begin(b));

        return a;
    }

    // Convert a ruint into a GMP integer
    template <size_t K>
    inline mpz_ptr ruint_to_mpz_t(mpz_ptr a, const ruint<K>& b) {
        // TODO Optimize...
        mpz_class r;
        RecInt::ruint_to_mpz(r, b);
        mpz_init_set(a, r.get_mpz_t()) ;
        return a;
    }

    // Convert a ruint into a GMP integer
    template <size_t K>
    inline ruint<K>& mpz_t_to_ruint(ruint<K>& a, mpz_srcptr b) {
        // TODO Optimize...
        mpz_class r(b);
        return mpz_to_ruint(a, r);
    }

    template <size_t K>
    inline ruint<K>::ruint(const char* b) {
        mpz_class m(b);
        mpz_to_ruint(*this, m);
    }

}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
