// =================================================================== //
// Givaro : Euler's phi function
//          Primitive roots.
//          RSA scheme.
// Time-stamp: <16 Apr 03 20:00:29 Jean-Guillaume.Dumas@imag.fr> 
// =================================================================== //

#ifndef _GIVARO_NUMTHEORY_
#define _GIVARO_NUMTHEORY_

#include <iostream>
#include "givaro/givinteger.h"
#include "givaro/givintprime.h"
#include "givaro/givintfactor.h"
#include "givaro/givrandom.h"

// =================================================================== //
// Givaro : Theoreme Chinois des Restes
// =================================================================== //
// const Integer tcr(Container<Integer>& Lr, Container<Integer>& Lm);
// const Integer& tcr(Integer& res, Integer& pm, Container<Integer>& Lr, Container<Integer>& Lm);

template<class RandIter = GivRandom>
class IntNumTheoDom : public IntFactorDom<RandIter> {
public:
    typedef typename IntFactorDom<RandIter>::Rep Rep;
    IntNumTheoDom(RandIter& g = *(new RandIter())) 
            :  IntFactorDom<RandIter>(g) {}
// =================================================================== //
// Euler's phi function
// =================================================================== //
    Rep& phi(Rep& r, const Rep& n) const ;
    template< template<class> class Container> Rep& phi(Rep& res, const Container<Rep>& Lf, const Rep& n) const ;


// =================================================================== //
// Primitive Root
// =================================================================== //
    Rep& prim_root(Rep&, const Rep&) const ;
    Rep& prim_root(Rep&, unsigned long&, const Rep&) const ;
    Rep& prim_root_of_prime(Rep&, const Rep&) const ;
    template<class Array> Rep& prim_root_of_prime(Rep& A, const Array& Lf, const Rep& phin, const Rep& n) const ;
    

    Rep& lowest_prim_root(Rep&, const Rep&) const ;
    bool is_prim_root(const Rep&, const Rep&) const ;
    Rep& order(Rep&, const Rep&, const Rep&) const ;
    bool isorder(const Rep&, const Rep&, const Rep&) const ;

// =================================================================== //
// Generalization of primitive roots for any modulus
// Primitive element, Primitive invertible
// Both functions coïncides except for m=8
// zeta : Order of a primitive element
// zeta_inv : Order of an invertible element
// Both functions coïncides except for m=8
// =================================================================== //
    Rep& prim_inv(Rep & , const Rep&) const ;
    Rep& prim_elem(Rep & , const Rep&) const ;
private:
    Rep& prim_base(Rep & , const Rep&) const ;
    Rep& zeta_base(Rep & , const Rep&) const ;
public:
    Rep& zeta_primpow(Rep & , const Rep&, unsigned long) const ;
    Rep& zeta_inv_primpow(Rep & , const Rep&, unsigned long) const ;
    Rep& zeta(Rep & , const Rep&) const ;
    Rep& zeta_inv(Rep & , const Rep&) const ;

// =================================================================== //
// Möbius function
// =================================================================== //

    template< template<class> class Container> short mobius(const Container<unsigned long>& lpow) const ;
    short mobius(const Rep& a) const;
};


#include "givaro/givintnumtheo.inl"

#endif
