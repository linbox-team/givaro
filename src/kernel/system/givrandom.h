// =================================================================== //
// Givaro : random generator
// a la Linbox ...
// Time-stamp: <18 Jan 07 15:04:09 Jean-Guillaume.Dumas@imag.fr> 
// =================================================================== //
#ifndef __GIVARO_RTANDOM__H__
#define __GIVARO_RTANDOM__H__

extern "C" {
# include <sys/time.h>
# include <sys/resource.h>
}

// -----------------------------------------------------
// Fishman, G.S. "Multiplicative congruential random
// number generators ..." Math. Comp. 54:331-344 (1990)

#if __GIVARO__DONOTUSE_longlong__
#define _GIVRAN_MULTIPLYER_ 950706376UL
#define _GIVRAN_MODULO_     2147483647UL
#else
#define _GIVRAN_MULTIPLYER_ 950706376ULL
#define _GIVRAN_MODULO_     2147483647ULL
#endif

class GivRandom {
    mutable unsigned long _seed;
public:
    typedef GivRandom random_generator;

    GivRandom(const unsigned long s = 0) 
            : _seed(s) {
        if (! s) {
		struct timeval tp;
		gettimeofday(&tp, 0) ;
		_seed = (long)(tp.tv_usec);
	}
        std::cerr << "seed : " << _seed << std::endl;
    }  
     
    GivRandom(const GivRandom& R) : _seed(R._seed) {}
    GivRandom& operator= (const GivRandom& R) { _seed = R._seed; return *this; }

    unsigned long seed() const { return _seed; }
    
    unsigned long operator() () const {
        return _seed = (unsigned long)( 
            (__GIVARO_INT64)_GIVRAN_MULTIPLYER_ 
            * (__GIVARO_INT64)_seed
            % (__GIVARO_INT64)_GIVRAN_MODULO_ );
    }

    template<class XXX> XXX& operator() (XXX& x) const {
        return x = (XXX)this->operator() ();
    }        
    
};

#endif
