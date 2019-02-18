// =============================================================
// Copyright(c)'1994-2013 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: J.-G. Dumas
// =============================================================

#ifndef __GIVARO_Spy_integer_H
#define __GIVARO_Spy_integer_H

#include "gmp++/gmp++_int.h"

namespace Givaro {
    /*! @internal
     * Spy structure to have access to protected members of Givaro::Integer.
     */
    struct SpyInteger {

        struct InHeritsInteger : public Integer {
        protected:
            friend struct SpyInteger;
        };

        static const InHeritsInteger::Rep* get_rep(const Integer& i) {
            return static_cast<const InHeritsInteger&>(i).get_rep();
        }

        static mpz_ptr get_mpz(Integer& i) {
            return static_cast<InHeritsInteger&>(i).get_mpz();
        }
        static mpz_ptr get_mpz(const Integer& i) {
            return const_cast<InHeritsInteger&>(static_cast<const InHeritsInteger&>(i)).get_mpz();
        }
        static mpz_srcptr get_mpz_const(const Integer& i) {
            return static_cast<const InHeritsInteger&>(i).get_mpz_const();
        }
    };
} // Givaro

#endif //__GIVARO_Spy_integer_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
