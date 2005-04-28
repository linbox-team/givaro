// ================================================================= //
// Givaro / Athapascan-1
// Irreducibily test
// Factorisations de  Polynomes dans Fp[X] :
//      Distinct Degree
//      Cantor-Zassenhaus
//      Berlekamp : in LinBox
// Time-stamp: <06 Jan 05 17:22:58 Jean-Guillaume.Dumas@imag.fr> 
// ================================================================= //
#ifndef _GIV_POLY1_FACTO_H_
#define _GIV_POLY1_FACTO_H_
#include <givaro/givrandom.h>
#include <givaro/givpoly1.h>

// template<class Domain, class StorageTag> class Poly1FactorDom {};

template<class Domain, class Tag, class RandIter = GivRandom>
class Poly1FactorDom : public Poly1Dom<Domain,Tag> {
protected:
    using Poly1Dom<Domain,Tag>::_domain;
    using Poly1Dom<Domain,Tag>::one;
    using Poly1Dom<Domain,Tag>::zero;
    typedef typename Poly1Dom<Domain,Tag>::Rep Rep;
    mutable RandIter _g;
public:
    typedef typename Poly1Dom<Domain,Tag>::element element;
    typedef RandIter random_generator;

        // Warning : there is a copy of the random Iterator ...
    Poly1FactorDom (Domain& d, const Indeter& X = Indeter(), const RandIter& g = RandIter() ) : Poly1Dom<Domain,Tag> (d,X), _g(g) {}
    Poly1FactorDom (const Poly1Dom<Domain,Tag>& P, const RandIter& g = RandIter()) : Poly1Dom<Domain,Tag> (P), _g(g) {}

// ---------------------------------------------------------------
// Splits a polynomial into prime factors of same degree
// ---------------------------------------------------------------

    template< template<class> class Container > void SplitFactor( 
        Container< Rep > & L
        , const Rep& G
        , Degree d
        , typename Domain::Residu_t MOD)  const ;
    
    template< template<class> class Container> void SplitFactor( 
        Container< Rep > & L
        , const Rep& G
        , Degree d) const {
        SplitFactor(L,G,d,_domain.residu());
    }


    Rep& SplitFactor(
        Rep& R
        , const Rep& G
        , Degree d
        , typename Domain::Residu_t MOD) const  ;
    
    Rep& SplitFactor(
        Rep& R
        , const Rep& G
        , Degree d) const  {
        return SplitFactor(R,G,d,_domain.residu() );
    }
    


// ---------------------------------------------------------------
// Splits a polynomial into divisors of homogenous prime factors
// ---------------------------------------------------------------

    template< template<class> class Container> void DistinctDegreeFactor(
        Container< Rep > & L
        , const Rep& f
        , typename Domain::Residu_t MOD)  const ;

    template< template<class> class Container> void DistinctDegreeFactor(
        Container< Rep > & L
        , const Rep& f)  const {
        DistinctDegreeFactor(L,f,_domain.residu());
    }            

// ---------------------------------------------------------------
// Cantor-Zassenhaus Polynomial factorization over Z/pZ
// ---------------------------------------------------------------

    template< template<class> class Container> void CZfactor( 
        Container< Rep > & Lf
        , Container< unsigned long > & Le
        , const Rep& f
        , typename Domain::Residu_t MOD)  const ;

    template< template<class> class Container> void CZfactor( 
        Container< Rep > & Lf
        , Container< unsigned long > & Le
        , const Rep& f )  const {
        CZfactor(Lf, Le, f,_domain.residu());
    }

// ---------------------------------------------------------------
// Gives one non-trivial factor of P if P is reducible
// returns P otherwise
// ---------------------------------------------------------------

    Rep& factor(
        Rep& W
        , const Rep& P
        , typename Domain::Residu_t MOD )  const ;
 
    Rep& factor(
        Rep& W
        , const Rep& P )  const {
        return factor(W,P,_domain.residu());
    }
    
            
// ---------------------------------------------------------------
// Irreducibility test
// ---------------------------------------------------------------

    bool is_irreducible(
        const Rep& P
        , typename Domain::Residu_t MOD )  const ;
    
    bool is_irreducible(const Rep& P )  const {
        return is_irreducible(P,_domain.residu());
    }
    
    bool is_irreducible2(
        const Rep& P
        , typename Domain::Residu_t MOD )  const ;
    
    bool is_irreducible2(const Rep& P )  const {
        return is_irreducible2(P,_domain.residu());
    }
    
            
// ---------------------------------------------------------------
// Irreducible polynomials
// ---------------------------------------------------------------
        /// random irreducible polynomial
    element& random_irreducible (element& P, Degree n) const ;
        /// random irreducible polynomial tries to be sparse
    element& creux_random_irreducible (element& P, Degree n) const ;

        /// random irreducible polynomial with X as primitive root
    element& ixe_irreducible (element& R, Degree n) const ;
        /// random irreducible polynomial with X as primitive root
    element& ixe_irreducible2 (element& R, Degree n) const ;

// ---------------------------------------------------------------
// Primitive polynomials
// ---------------------------------------------------------------

    IntegerDom::element order(const Rep& P, const Rep& F) const ;
    

    bool is_prim_root( const Rep& P, const Rep& F) const ;
    
    Rep& random_prim_root(Rep& P, Rep& R, Degree n) const ;
    

    Rep& give_random_prim_root(Rep& R, const Rep& F) const ;
    Rep& give_prim_root(Rep& R, const Rep& F) const ;

};



#include "givaro/givpoly1factor.inl" 
#include "givaro/givpoly1proot.inl" 
#endif 
