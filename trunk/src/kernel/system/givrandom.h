// =================================================================== //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Givaro : random generator
// a la Linbox ...
// Time-stamp: <13 Jul 07 14:40:27 Jean-Guillaume.Dumas@imag.fr>
// =================================================================== //

/*! @file system/givrandom.h
 * @ingroup system
 * @brief NO DOC
 * @bib Fishman, G.S. "Multiplicative congruential random
 * number generators ..." Math. Comp. 54:331-344 (1990)
 */

#ifndef __GIVARO_random_H
#define __GIVARO_random_H
#include <givaro/givconfig.h>

extern "C" {
# include <sys/time.h>
# include <sys/resource.h>
}


#if __GIVARO__DONOTUSE_longlong__
#define _GIVRAN_MULTIPLYER_ 950706376UL
#define _GIVRAN_MODULO_     2147483647UL
#else
#define _GIVRAN_MULTIPLYER_ 950706376ULL
#define _GIVRAN_MODULO_     2147483647ULL
#endif

namespace Givaro {

class GivRandom {
    mutable unsigned long _seed;
public:
    typedef GivRandom random_generator;

    GivRandom(const unsigned long s = 0)
            : _seed(s)
    {
        if (! s) {
		struct timeval tp;
		gettimeofday(&tp, 0) ;
		_seed = (unsigned long)(tp.tv_usec);
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

    unsigned long seed() const
    {
	    return _seed;
    }

// #if defined(__GIVARO_INT64)
#if 1
    unsigned long operator() () const
    {
        return _seed = (unsigned long)(
            (int64_t)_GIVRAN_MULTIPLYER_
            * (int64_t)_seed
            % (int64_t)_GIVRAN_MODULO_ );
    }
#else
    unsigned long operator() () const
    {
        return _seed = (unsigned long)(
            (unsigned long)_GIVRAN_MULTIPLYER_
            * _seed
            % (unsigned long)_GIVRAN_MODULO_ );
    }
#endif

    template<class XXX> XXX& operator() (XXX& x) const
    {
        return x = (XXX)this->operator() ();
    }

};

} // namespace Givaro

#endif // __GIVARO_random_H
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
