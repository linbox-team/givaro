// =================================================================== //
// Givaro : random generator
// a la Linbox ...
// Time-stamp: <26 Apr 01 18:15:28 Jean-Guillaume.Dumas@imag.fr> 
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
#define _GIVRAN_MULTIPLYER_ 950706376
#define _GIVRAN_MODULO_     2147483647

class GivRandom {
    unsigned long _seed;
public:
    typedef GivRandom random_generator;

    GivRandom(const unsigned long s = 0) 
            : _seed(s) {
        if (! s) {
		struct timeval tp;
		gettimeofday(&tp, 0) ;
		_seed = (long)(tp.tv_usec);
	}
    }  
     
    GivRandom(const GivRandom& R) : _seed(R._seed) {}
    GivRandom& operator= (const GivRandom& R) { _seed = R._seed; return *this; }

    unsigned long seed() { return _seed; }
    
    unsigned long operator() () {    
        return _seed = (unsigned long) ( 
            (unsigned long long)_GIVRAN_MULTIPLYER_ 
            * (unsigned long long)(_seed) 
            % (unsigned long long)_GIVRAN_MODULO_ );
    }

    template<class XXX> XXX& operator() (XXX& x) {
        return x = this->operator() ();
    }        
    
};

#endif
