// Copyright(c)'1994-2011 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// -------------------------------------------------------------
// Time-stamp: <12 Jun 15 16:48:44 Jean-Guillaume.Dumas@imag.fr>
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
#include <givaro/gfq.h>
#include <givaro/givpoly1.h>
#include <givaro/givinterpgeom.h>
#include <givaro/givinterpgeom-multip.h>


using namespace Givaro;

#ifdef __GIVARO_DEBUG
long long TTcount = 0;
#endif

struct BlackBoxPolynomial {
    typedef Poly1Dom<GFqDom<int64_t>,Dense > PolDom_t;
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


struct BlackBoxVectorOfPolynomial {
    typedef Poly1Dom<GFqDom<int64_t>,Dense > PolDom_t;
    typedef PolDom_t::Element	    Polynomial;
    typedef PolDom_t::Type_t	    Type_t;
    typedef std::vector< Type_t > Vect_t;
    typedef std::vector< Polynomial > VectPoly_t;

    const PolDom_t& _PD;
    VectPoly_t _func;


    BlackBoxVectorOfPolynomial(const PolDom_t& pold, VectPoly_t& P)
    : _PD(pold), _func(P)
    {}

    Vect_t& operator()(Vect_t& vi, const Type_t& xi) const {
        vi.resize(_func.size());
        Vect_t::iterator iter_vi = vi.begin();
        VectPoly_t::const_iterator iter_func = _func.begin();
        for( ; iter_func != _func.end(); ++iter_func, ++iter_vi)
            _PD.eval(*iter_vi, *iter_func, xi);
        return vi;
    }

};


template<class Interp, class RandIter>
bool TestOneField(Interp& FD, RandIter& generator, size_t degmax) {

    typename Interp::Polynomial_t nouv, prec;
    FD.getpoldomain().random(generator, prec, Degree((int64_t)degmax) );

#ifdef __GIVARO_DEBUG
    FD.getpoldomain().write(std::cout << "Precedent: ", prec) << std::endl;
#endif

    BlackBoxPolynomial TestBB(FD.getpoldomain(), prec);

    FD.initialize(TestBB);
    for(size_t i=1; i<=degmax; ++i)
        FD(TestBB);


    FD.interpolator(nouv);
#ifdef __GIVARO_DEBUG
    FD.getpoldomain().write(std::cout << "Nouveau: ", nouv) << std::endl;
#endif

    if (! FD.getpoldomain().areEqual(prec,nouv)) {
#ifdef __GIVARO_DEBUG
        std::cout << "ERREUR: Precedent != Nouveau"  << std::endl;
#endif
        return false;
    }

#ifdef __GIVARO_DEBUG
    ++TTcount;
#endif
    return true;
}


template<class Interp, class RandIter>
bool TestOneFieldVect(Interp& FD, RandIter& generator, size_t degmax, size_t numpol) {

    std::vector<typename Interp::Polynomial_t> nouv(numpol), prec(numpol);
    for(size_t i=0; i< numpol; ++i) {
        FD.getpoldomain().random(generator, prec[i], Degree((int64_t)degmax) );
#ifdef __GIVARO_DEBUG
        FD.getpoldomain().write(std::cout << "Precedent[" << i << "]: ", prec[i]) << std::endl;
#endif
    }

    BlackBoxVectorOfPolynomial TestBB(FD.getpoldomain(), prec);

    FD.initialize(TestBB);
    for(size_t i=1; i<=degmax; ++i)
        FD(TestBB);


    FD.interpolator(nouv);

#ifdef __GIVARO_DEBUG
    for(size_t i=0; i< numpol; ++i)
        FD.getpoldomain().write(std::cout << "Nouveau[" << i << "]: ", nouv[i]) << std::endl;
#endif



    for(size_t i=0; i< numpol; ++i)
        if (! FD.getpoldomain().areEqual(prec[i],nouv[i])) {
#ifdef __GIVARO_DEBUG
            std::cout << "ERREUR: Precedent[" << i << "] != Nouveau[" << i << ']'  << std::endl;
#endif
            return false;
        }

#ifdef __GIVARO_DEBUG
    ++TTcount;
#endif
    return true;
}


int main(int argc, char ** argv) {
    // argv[1] : modulo
    // argv[2] : deg max
    // argv[3] : num poly max
    // argv[4] : seed

    GFqDom<int64_t>::Residu_t MOD = (argc>1 ? (GFqDom<int64_t>::Residu_t) atoi(argv[1]) : 101);
    size_t degmax = (argc>2 ? (size_t)atoi(argv[2]) : 20);
    size_t numpolmax = (argc>3 ? (size_t)atoi(argv[3]) : 15);

    size_t seed = (argc>4?(size_t)atoi(argv[4]):(size_t)BaseTimer::seed());
#ifdef __GIVARO_DEBUG
    std::cerr << "seed: " << seed << std::endl;
#endif
    GivRandom generator(seed);

    GFqDom<int64_t> F(MOD);
    NewtonInterpGeom< GFqDom<int64_t> > FD(F, Indeter("X"));
    NewtonInterpGeomMultip< GFqDom<int64_t> > FDM(F, Indeter("X"));

    bool success = true;

    for(size_t d=1; d<degmax ; ++d)
        for(size_t i=0; i<50; ++i)
            success &= TestOneField(FD, generator, d);

    for(size_t p=1; p<numpolmax ; ++p)
        for(size_t d=1; d<degmax ; ++d)
            for(size_t i=0; i<5; ++i)
                success &= TestOneFieldVect(FDM, generator, d, p);

#ifdef __GIVARO_DEBUG
    if (! success) {
        std::cerr << "Error: " << seed << std::endl;
    } else {
        std::cerr << "Success:" << TTcount << std::endl;
    }
#endif

    return (! success);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
