// =================================================================== //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <10 Dec 21 13:22:09 Jean-Guillaume.Dumas@imag.fr>
// =================================================================== //

/*! @file givrandom.h
 * @ingroup system
 * @brief NO DOC
 * @bib
 *   - Fishman, GS <i>Multiplicative congruential random
 * number generators...</i> Math. Comp. 54:331-344 (1990).
 *
 */

#ifndef __GIVARO_random_H
#define __GIVARO_random_H
#include <givaro/givconfig.h>
#include <givaro/udl.h>
#include <givaro/givtimer.h>

extern "C" {
# include <sys/time.h>
# include <sys/resource.h>
}

#define _GIVRAN_MULTIPLYER_ 950706376_ui64
#define _GIVRAN_MODULO_     2147483647_ui64
#define _GIVRAN_HALFMOD_    1073741824_ui64

namespace Givaro {

    //! GivRandom
    class GivRandom {
        mutable uint64_t _seed;
    public:
        typedef GivRandom random_generator;
        typedef uint64_t random_t;
        typedef uint64_t seed_t;

        GivRandom(const uint64_t s = 0)
        : _seed(s)
        {
            while (! _seed) {
                _seed = (uint64_t)BaseTimer::seed();
            }
        }

        GivRandom(const GivRandom& R) :
            _seed(R._seed)
        {}

        GivRandom& operator= (const GivRandom& R)
        {
            _seed = R._seed;
            return *this;
        }

        uint64_t seed() const
        {
            return _seed;
        }

        uint64_t max_rand() const
        {
            return _GIVRAN_MODULO_;
        }

        uint64_t operator() () const
        {
            return _seed = (uint64_t)(
                                      (int64_t)_GIVRAN_MULTIPLYER_
                                      * (int64_t)_seed
                                      % (int64_t)_GIVRAN_MODULO_ );
        }

        bool brand() const
        {
            return !(this->operator()() & _GIVRAN_HALFMOD_);
        }

        template<class XXX> XXX& operator() (XXX& x) const
        {
            return x = (XXX)this->operator() ();
        }

    };

} // namespace Givaro

#endif // __GIVARO_random_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
