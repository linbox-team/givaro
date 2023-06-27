// ================================================================= //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <19 Apr 23 15:53:12 Jean-Guillaume.Dumas@imag.fr>
// ================================================================= //

/** @file givpoly1factor.h
 * @ingroup poly1
 * @brief NO DOC
 * - Irreducibily test
 * - Factorisations de  Polynomes dans Fp[X] :
 * - Distinct Degree
 * - Cantor-Zassenhaus
 * - Berlekamp: moved in LinBox
 * .
 */


#ifndef __GIVARO_poly1_facto_H
#define __GIVARO_poly1_facto_H
#include <givaro/givrandom.h>
#include <givaro/givpoly1.h>


namespace Givaro {

    // template<class Domain, class StorageTag> class Poly1FactorDom {};

    //! Poly1FactorDom
    template<class Domain, class Tag = Dense, class RandomIterator = GivRandom>
    class Poly1FactorDom : public Poly1Dom<Domain,Tag> {
    protected:
        using               Poly1Dom<Domain,Tag>::_domain;
        typedef typename Poly1Dom<Domain,Tag>::Rep    Rep;
        mutable          RandomIterator                      _g;
    public:
        using                                                 Poly1Dom<Domain,Tag>::one;
        using                                                Poly1Dom<Domain,Tag>::zero;
        typedef typename Poly1Dom<Domain,Tag>::Element                          Element;
        typedef          RandomIterator                                      random_generator;
        typedef typename Domain::Residu_t Residu_t;
        typedef typename Signed_Trait<typename Domain::Element>::signed_type   Element_t;
        // using typename Poly1Dom<Domain,Tag>::Type_t; //used for assign (clang complains)
        typedef typename Poly1Dom<Domain,Tag>::Type_t Type_t; //used for assign (clang complains)
        // typedef typename Domain::Residu_t                                      Residu_t;

        Poly1FactorDom () {}

        //! @warning  there is a copy of the random Iterator ...
        Poly1FactorDom (const Domain& d, const Indeter& X = Indeter(), const RandomIterator& g = RandomIterator() ) :
            Poly1Dom<Domain,Tag> (d,X), _g(g)
        {}

        Poly1FactorDom (const Poly1Dom<Domain,Tag>& P, const RandomIterator& g = RandomIterator()) :
            Poly1Dom<Domain,Tag> (P), _g(g)
        {}

        // ---------------------------------------------------------------
        // Splits a polynomial into prime factors of same degree
        // ---------------------------------------------------------------

        template< template<class, class> class Container, template<class> class Alloc >
        void SplitFactor( Container< Rep, Alloc<Rep> > & L
                          , const Rep& G
                          , Degree d
                          , Residu_t MOD)  const ;

        template< template<class, class> class Container, template <class> class Alloc>
        void SplitFactor( Container< Rep, Alloc<Rep> > & L
                          , const Rep& G
                          , Degree d) const {
            SplitFactor(L,G,d,_domain.residu());
        }


        Rep& SplitFactor(
                         Rep& R
                         , const Rep& G
                         , Degree d
                         , Residu_t MOD) const  ;

        Rep& SplitFactor(
                         Rep& R
                         , const Rep& G
                         , Degree d) const  {
            return SplitFactor(R,G,d,_domain.residu() );
        }



        // ---------------------------------------------------------------
        // Splits a polynomial into divisors of homogenous prime factors
        // ---------------------------------------------------------------

        template< template<class, class> class Container, template<class> class Alloc>
        void DistinctDegreeFactor(Container< Rep, Alloc<Rep> > & L
                                  , const Rep& f
                                  , Residu_t MOD)  const ;

        template< template<class, class> class Container, template <class> class Alloc>
        void DistinctDegreeFactor( Container< Rep, Alloc<Rep> > & L
                                   , const Rep& f)  const {
            DistinctDegreeFactor(L,f,_domain.residu());
        }

        // ---------------------------------------------------------------
        // Cantor-Zassenhaus Polynomial factorization over Z/pZ
        // ---------------------------------------------------------------

        template< template<class, class> class Container, template <class> class Alloc>
        void CZfactor( Container< Rep, Alloc<Rep> > & Lf,
                       Container< uint64_t, Alloc<uint64_t> > & Le,
                       const Rep& f,
                       Residu_t MOD)  const ;

        template< template<class, class> class Container, template <class> class Alloc>
        void CZfactor( Container< Rep, Alloc<Rep> > & Lf,
                       Container< uint64_t, Alloc<uint64_t> > & Le,
                       const Rep& f )  const {
            CZfactor(Lf, Le, f,_domain.residu());
        }

        // ---------------------------------------------------------------
        // Gives one non-trivial factor of P if P is reducible
        // returns P otherwise
        // ---------------------------------------------------------------

        Rep& factor(
                    Rep& W
                    , const Rep& P
                    , Residu_t MOD )  const ;

        Rep& factor(
                    Rep& W
                    , const Rep& P )  const {
            return factor(W,P,_domain.residu());
        }


        template<class PolyCont>
        PolyCont& factor (PolyCont& factors,
                          std::vector<uint64_t>& exp,
                          const Rep& P) {
            // JGD 02.03.2012 : to be refactored
            // at least without pointers ...
            CZfactor(factors, exp, P); // Cantor-Zassenhaus factorization
            return factors;
        }



        // ---------------------------------------------------------------
        // Irreducibility test
        // ---------------------------------------------------------------

        bool is_irreducible( const Rep& P
                             , Residu_t MOD )  const ;

        bool is_irreducible( const Rep& P )  const
        {
            return is_irreducible(P,_domain.residu());
        }

        bool is_irreducible2( const Rep& P
                              , Residu_t MOD )  const ;

        bool is_irreducible2( const Rep& P )  const
        {
            return is_irreducible2(P,_domain.residu());
        }


        // ---------------------------------------------------------------
        // Irreducible polynomials
        // ---------------------------------------------------------------
        /// random irreducible polynomial
        Element& random_irreducible (Element& P, Degree n) const ;
        /// random irreducible polynomial tries to be sparse
        Element& creux_random_irreducible (Element& P, Degree n) const ;

        /// random irreducible polynomial with X as primitive root
        Element& ixe_irreducible (Element& R, Degree n) const ;
        /// random irreducible polynomial with X as primitive root
        Element& ixe_irreducible2 (Element& R, Degree n) const ;

        // ---------------------------------------------------------------
        // Primitive polynomials
        // ---------------------------------------------------------------

        IntegerDom::Element order(const Rep& P, const Rep& F) const ;


        bool is_prim_root( const Rep& P, const Rep& F) const ;

        Rep& random_prim_root(Rep& P, Rep& R, Degree n) const ;


        Rep& give_random_prim_root(Rep& R, const Rep& F) const ;
        Rep& give_prim_root(Rep& R, const Rep& F) const ;

    protected :
        template<class Residue>
        bool find_irred_binomial(Element &R, Degree n, Residue MOD) const;
        bool find_irred_binomial(Element &R, Degree n, bool MOD) const;
        template<class Residue>
        bool find_irred_binomial(Element &R, Degree n, Residue MOD, Element IXE) const;
        bool find_irred_binomial(Element &R, Degree n, bool MOD, Element IXE) const;
        template<class Residue>
        bool find_irred_binomial2(Element &R, Degree n, Residue MOD, Element IXE) const;
        bool find_irred_binomial2(Element &R, Degree n, bool MOD, Element IXE) const;

        template<class Residue>
        bool find_irred_trinomial(Element &R, Degree n, Residue MOD) const;
        bool find_irred_trinomial(Element &R, Degree n, bool MOD) const;
        template<class Residue>
        bool find_irred_trinomial(Element &R, Degree n, Residue MOD, Element IXE) const;
        bool find_irred_trinomial(Element &R, Degree n, bool MOD, Element IXE) const;
        template<class Residue>
        bool find_irred_trinomial2(Element &R, Degree n, Residue MOD, Element IXE) const;
        bool find_irred_trinomial2(Element &R, Degree n, bool MOD, Element IXE) const;

        template<class Residue>
        bool find_irred_randomial(Element &R, Degree n, Residue MOD) const;
        bool find_irred_randomial(Element &R, Degree n, bool MOD) const;
        template<class Residue>
        bool find_irred_randomial(Element &R, Degree n, Residue MOD, Element IXE) const;
        bool find_irred_randomial(Element &R, Degree n, bool MOD, Element IXE) const;
        template<class Residue>
        bool find_irred_randomial2(Element &R, Degree n, Residue MOD, Element IXE) const;
        bool find_irred_randomial2(Element &R, Degree n, bool MOD, Element IXE) const;

    };

} // Givaro

#include "givaro/givpoly1factor.inl"
#include "givaro/givpoly1proot.inl"
#endif // __GIVARO_poly1_facto_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
