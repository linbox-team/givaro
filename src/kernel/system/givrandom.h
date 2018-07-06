// =================================================================== //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <28 Sep 16 12:03:51 Jean-Guillaume.Dumas@imag.fr>
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
#include <random>
#include <givaro/givconfig.h>
#include <givaro/udl.h>
#include <givaro/givtimer.h>

extern "C" {
# include <sys/time.h>
# include <sys/resource.h>
}

#define _GIVRAN_MULTIPLYER_ 950706376_ui64
#define _GIVRAN_MODULO_     2147483647_ui64

namespace Givaro {

	template <class EngineType>
    class GivRandomGeneric :public EngineType {
    public:
		typedef EngineType Engine;
        typedef GivRandomGeneric<Engine> random_generator;

		/** Initializes the PRNG with the given seed.
		 * A seed of zero means to get a true-random seed value.
		 */
        GivRandomGeneric(const uint64_t s = 0) :
			Engine(s)
		{
			if (!s) {
				std::random_device rd;
				Engine::seed(rd());
			}
		}

		uint64_t operator() ()
		{
			return Engine::operator()();
		}

        template<class XXX> XXX& operator() (XXX& x)
		{
			return x = (XXX)this->operator() ();
		}

    };

	using GivRandom = GivRandomGeneric<
		std::linear_congruential_engine
		<uint64_t, _GIVRAN_MULTIPLYER_, 0, _GIVRAN_MODULO_>>;

} // namespace Givaro

#endif // __GIVARO_random_H
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
