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
#include <gmpxx.h>
#include <givaro/givconfig.h>
#include <random>
#include <cstring>
#include <type_traits>
#include <givaro/udl.h>
#include <givaro/givtimer.h>

extern "C" {
# include <sys/time.h>
# include <sys/resource.h>
}

#define _GIVRAN_MULTIPLYER_ 950706376_ui64
#define _GIVRAN_MODULO_     2147483647_ui64

/* Define a macro that allows randomized functions to silently default to the
 * global RNG state, which should be deprecated in future code.
 */
#define _GIVARO_RAND_LEGACY // XXX FIXME should be a config option
#ifdef _GIVARO_RAND_LEGACY
#define _GIVRAN_DEFARG = default_givrand()
#else
#define _GIVRAN_DEFARG
#endif

namespace Givaro {

	// forward declaration
	class Integer;

	/** @brief Replacement for gmp_randclass that allows move constructor.
	 */
	class GivIntRand {
	protected:
		gmp_randstate_t _state;
		bool _initialized;

		// subclasses must initialize the state!
		GivIntRand() :_initialized(false) {}

	public:
		virtual ~GivIntRand()
		{
			if (_initialized) gmp_randclear(_state);
		}

		// disallow copying
		GivIntRand(const GivIntRand&) = delete;
		GivIntRand& operator= (const GivIntRand&) = delete;

		// allow moving
		GivIntRand(GivIntRand&& other) noexcept
		{
			std::memcpy(&_state, &other._state, sizeof _state);
			std::memset(&other._state, 0, sizeof _state);
			other._initialized = false;
		}
		GivIntRand& operator= (GivIntRand&& other) noexcept
		{
			if (_initialized) gmp_randclear(_state);
			std::memcpy(&_state, &other._state, sizeof(_state));
			std::memset(&other._state, 0, sizeof _state);
			other._initialized = false;
			return *this;
		}

		inline void seed(const uint64_t& s)
		{
			gmp_randseed_ui(_state, s);
		}

		// definition is in gmp++_int_rand.inl
		inline void seed(const Integer& s);

		// result <$- {0,...,m-1}
		// definition is in gmp++_int_rand.inl
		inline void get_z_range(Integer& result, const Integer& m);

		// result <$- {0,...,2^bits-1}
		// definition is in gmp++_int_rand.inl
		inline void get_z_bits(Integer& result, uint64_t bits);
	};

	/** @brief Linear congruential generator for Integers
	 */
	struct GivIntRandLC :public GivIntRand {
		static const mp_bitcnt_t _BITLEN = 64;

		template <typename SeedType=uint64_t>
		GivIntRandLC (const SeedType& seed=0)
		{
			gmp_randinit_lc_2exp_size(this->_state, _BITLEN);
			this->_initialized = true;
			this->seed(seed);
		}

		GivIntRandLC (GivIntRandLC&&) = default;
		GivIntRandLC (const GivIntRandLC&) = delete;
		GivIntRandLC& operator= (GivIntRandLC&&) = default;
		GivIntRandLC& operator= (const GivIntRandLC&) = delete;
	};

	/** @brief Mersenne Twister generator for Integers
	 */
	struct GivIntRandMT :public GivIntRand {
		template <typename SeedType=uint64_t>
		GivIntRandMT (const SeedType& seed=0)
		{
			gmp_randinit_mt(this->_state);
			this->_initialized = true;
			this->seed(seed);
		}

		GivIntRandMT (GivIntRandMT&&) = default;
		GivIntRandMT (const GivIntRandMT&) = delete;
		GivIntRandMT& operator= (GivIntRandMT&&) = default;
		GivIntRandMT& operator= (const GivIntRandMT&) = delete;
	};

	template <class EngineType, class IntRNGType>
    class GivRandomGeneric {
    public:
		typedef EngineType Engine;
		typedef IntRNGType IntegerRNG;
        typedef GivRandomGeneric<Engine,IntRNGType> random_generator;

	private:
		Engine _engine;
		IntegerRNG _intrng;

	public:
		/** Initializes the PRNG with the given seed.
		 * A seed of zero means to get a true-random seed value.
		 */
		template <typename SeedType=uint64_t>
        GivRandomGeneric(const SeedType& s = 0)
		{
			seed(s);
		}

#ifdef _GIVARO_RAND_LEGACY
		/** Copy constructor that creates a new PRNG seeded from the given one.
		 * Note that the copy will *not* produce the same sequence as the
		 * original.
		 * New code should explicitly use the newSeed() method instead.
		 */
		GivRandomGeneric(GivRandomGeneric& other)
		{
			seed(other.newSeed());
		}
#else
		GivRandomGeneric(GivRandomGeneric&) = delete;
#endif
		GivRandomGeneric(const GivRandomGeneric&) = delete;

		/** Move constructor. */
		GivRandomGeneric(GivRandomGeneric&&) = default;

		/** Move assignment. */
		GivRandomGeneric& operator= (GivRandomGeneric&&) = default;

		uint64_t operator() ()
		{
			return _engine();
		}

        template<class XXX> XXX& operator() (XXX& x)
		{
			return x = (XXX)_engine();
		}

		void seed(uint64_t s = 0) {
			while (!s) {
				std::random_device rd;
				s = rd();
			}
			_engine.seed(s);
			_intrng.seed(s);
		}

		// definition is in gmp++_int_rand.inl
		void seed(const Integer& s);

		/** Gets a new seed value suitable for instantiating a
		 * pseudo-independent RNG.
		 */
		uint64_t newSeed()
		{
			return this->operator()() ^ UINT64_C(0x9a52dc7e049f858d);
		}

		Engine& getEngine() { return _engine; }
		GivIntRand& getIntRand() { return _intrng; }
    };

	/** @brief Linear Congruential Generator RNG.
	 * This prioritizes speed over randomness quality.
	 */
	using GivRandomLC = GivRandomGeneric<
		std::linear_congruential_engine
			<uint64_t, _GIVRAN_MULTIPLYER_, 0, _GIVRAN_MODULO_>,
	    GivIntRandLC>;

	/** @brief Mersenne Twister RNG.
	 * This still has very good performance, but with better randomness
	 * properties, compared to LCG.
	 */
	using GivRandomMT = GivRandomGeneric<
		std::mt19937_64, GivIntRandMT>;

	/** @brief RNG using default settings. */
	using GivRandom = GivRandomLC;

	inline GivRandom& default_givrand() {
		static GivRandom gr;
		return gr;
	}

	/** @brief Type trait for any class following the GivRandom design.
	 */
	template <class T>
	struct IsGivRand :std::false_type {};

	template <class EngineType, class IntRNGType>
	struct IsGivRand<GivRandomGeneric<EngineType,IntRNGType>> :std::true_type {};

} // namespace Givaro

#endif // __GIVARO_random_H
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
