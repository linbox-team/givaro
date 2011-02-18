// Copyright(c)'1994-2011 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software. 
// see the COPYRIGHT file for more details.
// -------------------------------------------------------------
// Time-stamp: <18 Feb 11 16:11:21 Jean-Guillaume.Dumas@imag.fr>
//                                                              
// Interpolation at geometric points                            
// see: Polynomial evaluation and interpolation on special sets 
//      of points,                                              
//      A. Bostan and E. Schost,                                
//      Journal of Complexity 21(4): 420-446, 2005.             
// -------------------------------------------------------------

#include <iostream>
#include <givaro/givrandom.h>
#include <givaro/givtimer.h>
#include <givaro/givgfq.h>
#include <givaro/givpoly1.h>
#include <givaro/givinterpgeom.h>


#ifdef GIVARO_DEBUG
long long TTcount = 0;
#endif

struct BlackBoxPolynomial {
    typedef Poly1Dom<GFqDom<long>,Dense > PolDom_t;
    typedef PolDom_t::Element	    Polynomial;
    typedef PolDom_t::Type_t	    Type_t;
    
    const PolDom_t& _PD;
    Polynomial _func;
    

    BlackBoxPolynomial(const PolDom_t& pold, Polynomial& P) 
            : _PD(pold), _func(P) 
        {}
    
    Type_t& operator()(Type_t& vi, const Type_t& xi) const {
        return _PD.eval(vi, _func, xi);
    }
    
};

template<class Interp, class RandIter>
bool TestOneField(Interp& FD, RandIter& generator, size_t degmax) {

    typename Interp::Polynomial_t nouv, prec;
    FD.getpoldomain().random(generator, prec, Degree(degmax) );

#ifdef GIVARO_DEBUG
    FD.getpoldomain().write(std::cout << "Precedent: ", prec) << std::endl;
#endif

    BlackBoxPolynomial TestBB(FD.getpoldomain(), prec);
    
    FD.initialize(TestBB);
    for(size_t i=1; i<=degmax; ++i)
        FD(TestBB);
        

    FD.interpolator(nouv);
#ifdef GIVARO_DEBUG
    FD.getpoldomain().write(std::cout << "Nouveau: ", nouv) << std::endl;
#endif 
    
    if (! FD.getpoldomain().areEqual(prec,nouv)) {
#ifdef GIVARO_DEBUG
        std::cout << "ERREUR: Precedent != Nouveau"  << std::endl;
#endif
        return false;
    }

#ifdef GIVARO_DEBUG
    ++TTcount;
#endif
    return true;
}


int main(int argc, char ** argv) {
        // argv[1] : modulo
        // argv[2] : deg max
        // argv[3] : seed

    GFqDom<long>::Residu_t MOD = (argc>1 ? atoi(argv[1]) : 101);
    size_t degmax = (argc>2 ? atoi(argv[2]) : 20);

    size_t seed = (argc>3?atoi(argv[3]):BaseTimer::seed());
#ifdef GIVARO_DEBUG
    std::cerr << "seed: " << seed << std::endl;
#endif
    GivRandom generator(seed);

    GFqDom<long> F(MOD);
    NewtonInterpGeom< GFqDom<long> > FD(F, Indeter("X"));

    bool success = true;
    for(size_t d=1; d<degmax ; ++d)
        for(size_t i=0; i<50; ++i)
            success &= TestOneField(FD, generator, d);

#ifdef GIVARO_DEBUG
    if (! success) {
        std::cerr << "Error: " << seed << std::endl;
    } else {
        std::cerr << "Success:" << TTcount << std::endl;
    }
#endif

    return (! success);
}
